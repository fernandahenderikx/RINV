%% process IFCB data Integrated hourly to improve SNR of rare particles
clearvars -except mdir procdata datadir figsdir
mdir=cd;
cd ..
cd data
datadir = cd;
cd ..\proc
procdata = cd;
cd(mdir)

load([datadir '\all_features_small.mat']);
load([procdata '\refdiam.mat']);

ifcb = data; clear data
ifcbTime    = datenum(ifcb.Date);
ifcbVolAn   = ifcb.Volume_analyzed.*10^-3; % cm3 to L
ifcbBiovol  = ifcb.Biovolume.*10^-9; % um3 to uL 
ifcbDiam    = ifcb.ESD; %um, from BIOVOLUME
%ifcbDiam    = ((ifcb.CSA./pi).^(1/2)).*2; %um, from Cross-sectional area
ifcbVolConc = ifcbBiovol./ifcbVolAn; % ul/L --> Accounts for Volume Analyzed

% get size classes - make it match LISST
% for each time step, get the # particles in each size class.
% size bins of LISST:

% Bubble check 
indBad1 = [23912:26651 252914:255348];
% bead check
indBad2 = find(ifcbTime<=datenum(2017,06,27,17,0,0))';
indBad3 = find(ifcbDiam==0)';
indBad = [indBad1 indBad2 indBad3];
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
binwidth  = bins(:,2) - bins(:,1);

% GET IFCB #/L and ul/L  at same size bins as LISST
[it,uu]        = unique(ifcbTime);

% GET IFCB #/L and ul/L, bin to every hour
temp           = datevec(min(ifcbTime));
startTime      = datenum([temp(1:4),0,0]);
temp           = datevec(max(ifcbTime));
endTime        = datenum([temp(1:3),temp(4)+1,0,0]);
timeHour       = [startTime:(1/24):endTime];

ifcbVolconcAggBinHr = NaN(length(timeHour)-1,1);
ifcbMedDiamBinHr  = NaN(length(timeHour)-1,1);
ifcbvolAnHr       = NaN(length(timeHour)-1,1);
ifcbPsdBinHr      = NaN(length(timeHour)-1,38);
ifcbVolconcBinHr  = NaN(length(timeHour)-1,38);
ifcbMeanDiameterBinHr = NaN(length(timeHour)-1,38);
for iHour = 1:length(timeHour)-1
    indTime = find(ifcbTime >= timeHour(iHour) & ifcbTime <= timeHour(iHour+1));
    ifcbvolAnHr(iHour)         = sum(unique(ifcbVolAn(indTime)));
    ifcbVolconcAggBinHr(iHour) = nansum(ifcbBiovol(indTime))/ifcbvolAnHr(iHour);
    
    datavolc                  = ifcbBiovol(indTime);
    datadiam                  = ifcbDiam(indTime);
    xx = find(datadiam>=3.98); % only calculate diameter for sizes > 4um
    ifcbMedDiamBinHr(iHour)   = nanmedian(datadiam(xx));   

    for iBin = 1:38 % good data start at 4.3um
        [ind,indd] = find((datadiam >= bins(iBin,1)) & (datadiam < bins(iBin,2)));
        if isempty(ind);continue;end
        ifcbPsdBinHr(iHour,iBin)     = length(ind)./ifcbvolAnHr(iHour); % #/L, account for vol analyzed
        ifcbVolconcBinHr(iHour,iBin) = nansum(datavolc(ind))/ifcbvolAnHr(iHour); % uL/L
        ifcbMeanDiameterBinHr(iHour,iBin)   = nanmean(datadiam(ind)); % this is the mean diameter of each bin to be used instead of geometric mean from bin size!

    end
end
ifcbTimeHr = timeHour(1:end-1);

% crazy values:
bad = [12 20 21 70 60 71 72 110:113 167 171:180 183 197 250:252 291:297 305 309 313 317 411:414];
ifcbPsdnBinHr(bad,:)=NaN;
ifcbVolconcBinHr(bad,:)=NaN;
ifcbPsdBinHr(bad,:)=NaN;
ifcbMedDiamBinHr(bad,:)=NaN;
ifcbVolconcAggBinHr(bad) = NaN;

ifcbPsdnBinHr       = ifcbPsdBinHr./repmat(binwidth',length(timeHour)-1,1);
ifcbVolconcBinHrCut = nansum(ifcbVolconcBinHr(:,refdiam:end),2); % uL/L for >4um
ifcbVolconcBinHrCut(ifcbVolconcBinHrCut==0)=NaN;

cd(procdata)
save ifcb_IntegratedHourly.mat ifcb*Hr* bins binwidth
cd(mdir)