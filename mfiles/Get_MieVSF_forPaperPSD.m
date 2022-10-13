%% get Mie VSF for LISST angles and IFCB PSD + extrapolated portion
close all;
clearvars -except procdata mdir datadir figsdir
mdir = cd;
cd ../proc
load('ifcb_IntegratedHourly.mat','binwidth','ifcbPsdnBinHr','ifcbTimeHr')
load('refdiam.mat','dias'); % dias var was obtained from weighted averaging and NOT geometric mean. Use this instead of bins(:,3) from above mat file for consistency with the paper.


%%
data    = ifcbPsdnBinHr;
gmt     = ifcbTimeHr;
slope   = 4.57;
refdiam = 21;
No      = 5.310886e+04;

% build full PSD:
% first, clear bins 1:20 for clarity (we don't trust IFCB data at those
% bins)
data(:,(1:refdiam - 1)) = NaN;

% assume 1.25-6um size fraction does not change over time, so assume a fixed slope and No
N02_6 = No.*((dias(1:refdiam - 1)./dias(refdiam)).^-slope);

% add those extrapolated portions to the data to make a full psd
fullpsdn = data;
fullpsdn(:,(1:refdiam - 1)) = repmat(N02_6,length(fullpsdn(:,1)),1);

% add same gaps as measured data jsut in case this biases things later on
ind                           = isnan(fullpsdn(:,refdiam));
fullpsdn(ind,(1:refdiam - 1)) = NaN;
fullpsd                       = fullpsdn.*binwidth';
clear data fullpsdn

% calcmie below requires uL/L
class_vol = (1/6)*pi*(dias.^3); %um^3, equivalent 4/3PIr^3 and 1/6PId^3
class_vol = class_vol*1e-9; %ul, 10^9 um3 in a ul
vdc       = fullpsd.*(repmat(class_vol,size(fullpsd,1),1));  % uL/L

% get rid of NaNs for now
ind             = find(isnan(vdc(:,1)));
vdc(ind,:)      = [];
vdc(isnan(vdc)) = 0;
gmt(ind)        = [];
fullpsd(ind,:)  = [];

%% Get M11 Mueller Matrix for a range of n, calculate VSF for LISST angles

% Set up inputs
criRIn  = (1.015:0.0001:1.0333); % real refractive indices - just doing this limited range here to save time because I saw where the results were converging...
criI    = 0.0003; % immaginary component of the refractive index - same as paper.
volDist = vdc;
wave    = 670; % wavelength of the LISST instrument in nm.
nmedia  = 1; % water refractive index. use 1 for relative cri

% Compute Volumes, get diams and binwidths in the correct units (meter)
diams     = dias*10^-6; % m
binWidth  = binwidth'*10^-6; % m
binVol    = ((1/6)*pi*(diams.^3)).*1e18.*1e-9; % m3 > um3 > uL

% Now convert volDist to integrative number distribution
nDist = volDist./binVol./binWidth.*1e3; % uL/L > #/L > #/L/m > # m-4


% Scattering angles
% Basically make sure at least 20 points in each ring angle, then put on
% top of the original 821 angles.
baseAng     = logspace(0,log10(200),33)'*0.1; % B100X:0.1, C100X:0.05
binAng(:,1) = baseAng(1:32);% The lower limits are the first 32 (of 33)
binAng(:,2) = baseAng(2:33);% The upper limits are the last 32 (of 33)
binAng(:,3) = sqrt(binAng(:,1).*binAng(:,2));% Midpoints.
binAng      = binAng./1.33; % in-air to in-water (assuming nmed is 1)
vsfAngLisst = [];
for iAng = 1:32
    vsfAngLisst = [vsfAngLisst,logspace(log10(baseAng(iAng)/1.33),log10(baseAng(iAng+1)/1.33),20)];
end

%Set up scattering angles and size parameter
vsfAng821  = [0:0.01:2,2.025:0.025:5,5.5:0.5:9.5,10:170,170.5:0.5:175,...
              175.025:0.025:178,178.01:0.01:180]';% Mie VSF angles
vsfAngFull = sort([vsfAngLisst(:);vsfAng821(:)]);
lambda     = wave*1E-9;
r          = diams./2;           % convert to radius
k          = 2*pi/lambda;     % wavenumber in medium nm
x          = k*r;                % size parameter

% get VSF
vsf_ifcb_ang14 = nan(length(criRIn),size(fullpsd,1));
vsf_ifcb_ang18 = nan(length(criRIn),size(fullpsd,1));
vsf_ifcb_ang23 = nan(length(criRIn),size(fullpsd,1));
vsf_ifcb_ang22 = nan(length(criRIn),size(fullpsd,1));
vsf_ifcb_ang21 = nan(length(criRIn),size(fullpsd,1));
vsf_ifcb_ang20 = nan(length(criRIn),size(fullpsd,1));
vsf_ifcb_ang24 = nan(length(criRIn),size(fullpsd,1));
vsf_ifcb_ang25 = nan(length(criRIn),size(fullpsd,1));
vsf_ifcb_ang26 = nan(length(criRIn),size(fullpsd,1));

tic
for iCri = 1:length(criRIn)
    cri    = complex(criRIn(iCri),criI);
    Qext = NaN(length(diams),1);
    Qsca = NaN(length(diams),1);
    q11  = NaN(length(diams),length(vsfAngFull));
    for iDiam = 1:length(x) % hmm..this is not vectorized, so relatively slow
        [~,~,QscaIn,QextIn,~,~,Q11In] = fastmie_mod(x(iDiam),cri,[],vsfAngFull*pi/180);
        Qext(iDiam) = QextIn;
        Qsca(iDiam) = QscaIn;
        q11(iDiam,:) = Q11In;
    end

    %cp     = (pi/4).*(trapz(diams,diams.^2.*Qext'.*nDist,2));
    %bp     = (pi/4).*(trapz(diams,diams.^2.*Qsca'.*nDist,2));

    % Get VSF (M11)
    y   = pi*diams'*nmedia/lambda;% Convert diameter to size parameter for proper integration
    c11 = 4*q11./(y.^2); % c1 = directional efficiencies of the scattering cross section

    m11Full = NaN(size(vsfAngFull,1),size(nDist,1));
    for iDist = 1:size(nDist,1) % Too lazy, just loop through
        m11Full(:,iDist)  = (pi/4).*trapz(diams,diams.^2.*c11'.*nDist(iDist,:),2); % C11 is multiplied by the geometrical cross-section (with diameter converted to size parameter space since that's how fastmie does it). Result is M11.
    end

    % LISST ring averaged VSF (P11, the normalized mueller scattering matrix)
    p11Lisst = NaN(32,size(nDist,1));
    for iAng = 1:32
        indAng = (vsfAngFull>=binAng(iAng,1)) & (vsfAngFull<=binAng(iAng,2));
        p11Lisst(iAng,:) = trapz(vsfAngFull(indAng),m11Full(indAng,:).*sind(vsfAngFull(indAng)))./...
            trapz(vsfAngFull(indAng),sind(vsfAngFull(indAng)));
    end
    p11angLisst = binAng(:,3); % 

    vsf_ifcb_ang14(iCri,:) = p11Lisst(14,:);
    vsf_ifcb_ang18(iCri,:) = p11Lisst(18,:);
    vsf_ifcb_ang23(iCri,:) = p11Lisst(23,:);
    vsf_ifcb_ang22(iCri,:) = p11Lisst(22,:);
    vsf_ifcb_ang21(iCri,:) = p11Lisst(21,:);
    vsf_ifcb_ang20(iCri,:) = p11Lisst(20,:);
    vsf_ifcb_ang24(iCri,:) = p11Lisst(24,:);
    vsf_ifcb_ang25(iCri,:) = p11Lisst(25,:);
    vsf_ifcb_ang26(iCri,:) = p11Lisst(26,:);
end
toc

nint = criRIn;
ang_target = binAng(:,3); % Should this be the third column instead?
save vsf_ifcb_Full_ALLn.mat vsf_* vdc gmt nint fullpsd ang_target slope No