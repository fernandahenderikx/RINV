%% 

mdir=cd;
cd ..
cd data
datadir = cd;
cd ..\proc
procdata = cd;
cd(mdir)

load([datadir '\IFCB_Biovolume_raw.mat']);
ifcbTime    = ifcb(:,1);
ifcbVolAn   = ifcb(:,4).*10^-3; % cm3 to L
ifcbBiovol  = ifcb(:,2).*10^-9; % um3 to uL
ifcbDiam    = ifcb(:,3); %um
ifcbVolConc = ifcbBiovol./ifcbVolAn; % ul/L --> Accounts for Volume Analyzed
% get size classes - make it match LISST
% for each time step, get the # particles in each size class.
% size bins of LISST:

% Bubble check 
indBad = [23912:26651 252914:255348];
%BAD = find(ifcbBiovol>0.0003);
ifcbTime(indBad)    = [];
ifcbVolAn(indBad)   = [];
ifcbBiovol(indBad)  =  [];
ifcbDiam(indBad)    =  [];
ifcbVolConc(indBad) =  []; 

% Get diameters (select for spheres or randomly shaped particles)
x         = 0.2023; % lower limit 
rho       = 200^(1/32);
bins(:,1) = x*rho.^([0:37]); % lower limit
bins(:,2) = x*rho.^([1:38]); % upper limit
bins(:,3) = sqrt(bins(:,1).*bins(:,2)); % mid-point
dias      = bins(:,3)';
class_vol = (1/6)*pi*(dias.^3); % um^3, equivalent 4/3PIr^3 and 1/6PId^3
class_vol = class_vol*1e-9; % ul, 10^9 um3 in a uL
binwidth  = bins(:,2) - bins(:,1);

% GET IFCB #/L and ul/L  at same size bins as LISST
[it,uu]        = unique(ifcbTime);

% We are interested in a cruise-averaged PSD fit, so Just bin everything
% together to improve SNR
ifcbPsdFull     = NaN(1,38);
ifcbCountFull   = NaN(1,38);
ifcbVolconcFull = NaN(1,38);
volanalyzedFull = sum(ifcbVolAn(uu));
for iBin = 1:38 % good data start at 4.3um
    [ind,indd] = find((ifcbDiam >= bins(iBin,1)) & (ifcbDiam < bins(iBin,2)));
    ifcbPsdFull(1,iBin)           = length(ind)./volanalyzedFull; % #/L, account for vol analyzed
    ifcbCountFull(1,iBin)         = length(ind); % #/L,get counts for stats
    ifcbVolconcFull(:,iBin) = nansum(ifcbVolConc(indd)); % uL/L
end
ifcbPsdnFull = ifcbPsdFull./binwidth';

% Get log-log regression of full cruise PSD
indBin   = [19:30]; % 19 is the first good bin; bin 30 is the last bin with at least ~10 samples...bins>30 are likely still underestimating counts at those sizes
mdl      = fitlm(log10(dias(indBin)')-log10(4.3279),log10(ifcbPsdnFull(indBin)'));
modelfun = @(b,x)b(1).*(x./4.3279).^b(2);
x1       = logspace(log10(0.2),log10(4.3279),50);
y1       = modelfun([10^mdl.Coefficients{1,1},mdl.Coefficients{2,1}],x1);
plot(x1,y1,'k-')
x1       = logspace(log10(3),log10(70),50);
y1       = modelfun([10^mdl.Coefficients{1,1},mdl.Coefficients{2,1}],x1);
figure;plot(log10(dias(indBin)'),log10(ifcbPsdnFull(indBin)'),'.')
hold on;plot(log10(x1),log10(y1),'.')

% Print out report
out = mdl.Coefficients{1:2,1:2};
fprintf('\nSlope: %0.4f, SE: %0.4f, STD: %0.4f\nN0:10^%0.4f, SE of exponent: %0.4f, STD of exponent: %0.4f\n',...
    -1*out(2),out(4),out(4)*sqrt(length(indBin)),out(1),out(3),out(3)*sqrt(length(indBin)))
fprintf('\nN0(+1sigma): %i\nN0: %i\nN0(-1sigma): %i\n\n',10^(out(1)+out(3)*sqrt(length(indBin))),...
    10^out(1),10^(out(1)-out(3)*sqrt(length(indBin))))

cd(procdata)
fileID = fopen('ifcbPSDfitReport.txt','w');
fprintf(fileID,'\nSlope: %0.4f, SE: %0.4f, STD: %0.4f\nN0:10^%0.4f, SE of exponent: %0.4f, STD of exponent: %0.4f\n',...
    -1*out(2),out(4),out(4)*sqrt(length(indBin)),out(1),out(3),out(3)*sqrt(length(indBin)))
fprintf(fileID,'\nN0(+1sigma): %i\nN0: %i\nN0(-1sigma): %i\n\n',10^(out(1)+out(3)*sqrt(length(indBin))),...
    10^out(1),10^(out(1)-out(3)*sqrt(length(indBin))))
fclose(fileID);

cd(mdir)