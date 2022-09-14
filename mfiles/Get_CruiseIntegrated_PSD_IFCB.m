%% 
mdir=cd;
cd ..
cd data
datadir = cd;
cd ..\proc
procdata = cd;
cd(mdir)

load([datadir '\all_features_small.mat']);
ifcb = data; clear data
ifcbTime    = datenum(ifcb.Date);
ifcbVolAn   = ifcb.Volume_analyzed.*10^-3; % cm3 to L
ifcbBiovol  = ifcb.Biovolume.*10^-9; % um3 to uL 
ifcbDiam    = ifcb.ESD; %um, from BIOVOLUME
%ifcbDiam    = ((ifcb.CSA./pi).^(1/2)).*2; %um, from Cross-sectional area

% get size classes - make it match LISST
% for each time step, get the # particles in each size class.
% size bins of LISST:

% Bubble check 
indBad1 = [23912:26651 252914:255348];
indBad2 = find(ifcbTime<=datenum(2017,06,27,17,0,0))';
% bead check
indBad3 = find(ifcbDiam==0)';
indBad = [indBad1 indBad2 indBad3];
%BAD = find(ifcbBiovol>0.0003);
ifcbTime(indBad)    = [];
ifcbVolAn(indBad)   = [];
ifcbBiovol(indBad)  = [];
ifcbDiam(indBad)    = [];
ifcb(indBad,:)=[];

% Get bin sizes 
x         = 0.2023; % lower limit 
rho       = 200^(1/32);
bins(:,1) = x*rho.^([0:37]); % lower limit
bins(:,2) = x*rho.^([1:38]); % upper limit
binwidth  = bins(:,2) - bins(:,1);



%% BIN PSD DATA, GET nominal mean diameter per bin, fit a power law with weights
[it,uu]        = unique(ifcbTime);

% We are interested in a cruise-averaged PSD fit, so Just bin everything
% together to improve SNR
ifcbPsdFull     = NaN(1,38);
ifcbCountFull   = NaN(1,38);
ifcbMeanDiameterBinFull = NaN(1,38);
volanalyzedFull = sum(ifcbVolAn(uu));
for iBin = 1:38 % good data start at 4.3um
    [ind,indd] = find((ifcbDiam >= bins(iBin,1)) & (ifcbDiam < bins(iBin,2)));
    ifcbPsdFull(1,iBin)                     = length(ind)./volanalyzedFull; % #/L, account for vol analyzed
    ifcbCountFull(1,iBin)                   = length(ind); % # per bin ,get counts for stats
    ifcbMeanDiameterBinFull(:,iBin)         = nanmean(ifcbDiam(ind)); % this is the mean diameter of each bin to be used instead of geometric mean from bin size!
end
ifcbPsdnFull = ifcbPsdFull./binwidth'; % hmmm this is weird, no? I use the full binwidth, but the mean diameter is the weighted average of the points. Should not be a problem if there are enough datapoints.


figure
loglog(ifcbMeanDiameterBinFull,ifcbPsdnFull,'.')

% first good IFCB diameter:
[ind, refdiam]= nanmax(ifcbPsdnFull);


refdiam_value = ifcbMeanDiameterBinFull(refdiam);
% Get log-log regression of full cruise PSD - assume weights ~ number of particles per bin
indBin   = [refdiam:30]; % refdiam is the first good bin; bin 30 is the last bin with at least ~10 samples...bins>30 are likely still underestimating counts at those sizes
mdl      = fitlm(log10(ifcbMeanDiameterBinFull(indBin)')-log10(refdiam_value),log10(ifcbPsdnFull(indBin)'),'Weights',ifcbCountFull(indBin));
modelfun = @(b,x)b(1).*(x./refdiam_value).^b(2);
x1       = logspace(log10(0.2),log10(refdiam_value),50);
y1       = modelfun([10^mdl.Coefficients{1,1},mdl.Coefficients{2,1}],x1);
plot(x1,y1,'k-')
x1       = logspace(log10(3),log10(70),50);
y1       = modelfun([10^mdl.Coefficients{1,1},mdl.Coefficients{2,1}],x1);
figure;plot(log10(ifcbMeanDiameterBinFull(indBin)'),log10(ifcbPsdnFull(indBin)'),'o','MarkerFaceColor','b')
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

save ifcb_cruiseIntegrated.mat 

cd(mdir)

%% get sphericity proxy - check that particles are generally "spherical"

r = ifcb.MajorAxis./ifcb.MinorAxis;
ecce = ifcb.MinorAxis./ifcb.MajorAxis;

figure
hist(ecce,1000)

figure
plot(ifcb.ESD,ecce,'.')

bins(:,3)=sqrt(bins(:,1).*bins(:,2)); % use geometric mean here just to bin more easily

idx=[];
for i = 1:length(bins(:,1)) % just to guarantee I'm picking up everything, and that way size of idx is == size of r
    if i ==38
    ind = find(ifcb.ESD>=bins(i,3));
    else
       ind = find(ifcb.ESD>=bins(i,3)&ifcb.ESD<bins(i+1,3));
    end
    idx = [idx; ones(length(ind),1).*bins(i,3)];
end

figure
boxplot(ecce,idx)
ylim([0 1])

close all