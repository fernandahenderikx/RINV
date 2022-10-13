%% get Mie VSF for LISST angles and IFCB PSD + extrapolated portion

% this takes a long time to run as is.

mdir=cd;
cd ../proc
load('ifcb_IntegratedHourly.mat','bins','binwidth','ifcbPsdnBinHr','ifcbTimeHr')
load('refdiam','dias')

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



%% all of the below jusst to get the q11 matrix

% Set up inputs
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


cri = complex(1.02,0.0003);
% pick one time point: sunrise of July 6 2017
j=185;

   
    q11  = NaN(length(diams),length(vsfAngFull));
    for iDiam = 1:length(x) % hmm..this is not vectorized, so relatively slow
        [~,~,~,~,~,~,Q11In] = fastmie_mod(x(iDiam),cri,[],vsfAngFull*pi/180);

        q11(iDiam,:) = Q11In;
    end

    % Get VSF
    y   = pi*diams'*nmedia/lambda;% Convert diameter to size parameter for proper integration
    c11 = 4*q11./(y.^2); % q11 efficiency accounting for size parameter

% look at the VSD per ring
x = diams.^2.*c11'.*nDist(j,:);
xx = x(842,:); % this is angle that corresponds to LISST angle #26
xxx = xx.*binwidth'; % prepare for "integration" below since x axis is not evenly spaced in linear space

figure
plot(diams.*10^6,cumsum(xxx)./sum(xxx),'.') % 85% of signal comes from particle sizes between ~2 and 14 micron
xlabel('diameter'); ylabel('normalized cumulative sum of contribution of different size classes to VSF ring 26')
%1.87806395230894	0.112687709016137
%13.696090759242	0.958276387739291




