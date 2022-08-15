% Load data from Fernanda
mdir=cd;
cd ..
cd data
datadir = cd;
cd ..\proc
procdata = cd;
cd(mdir)


load([datadir '\IFCB_Biovolume_raw.mat']);
ifcbTime    = ifcb(:,1);
clear ifcb
load([datadir '\eco_proc.mat'],'eco_corr','time_min');
load([datadir '\proc_beamc_NEW.mat'],'beamc_corr','min_time');
bbp     = eco_corr(:,2);
bbpTime = time_min';
cp      = beamc_corr';
cpTime  = min_time';
clear eco_corr time_min beamc_corr min_time

% Take out NaNs already existing
indBad = isnan(bbp);
bbp(indBad) = [];
bbpTime(indBad) = [];
indBad = isnan(cp);
cp(indBad) = [];
cpTime(indBad) = [];

% Initial Figure
figsiz = [0 0 10 6].*1.6;
figure('Units','inches','Position',figsiz,...
    'PaperSize',figsiz(3:4),'PaperPosition',figsiz);
% tlo = tiledlayout(2,1,'TileSpacing','none','Padding','compact');
% nexttile;
yyaxis left
plot(cpTime,cp,'.'); hold on
%nexttile;
yyaxis right
plot(bbpTime,bbp,'.')

% Check 4 sigma outliers
indBad = 1;
while sum(indBad) > 0
indBad = (cp < nanmean(cp) - 4*nanstd(cp)) | (cp > nanmean(cp) + 4*nanstd(cp));
cp(indBad) = NaN;
end
indBad = 1;
while sum(indBad) > 0
indBad = (bbp < nanmean(bbp) - 4*nanstd(bbp)) | (bbp > nanmean(bbp) + 4*nanstd(bbp));
bbp(indBad) = NaN;
end

% Report out
rateGood = [sum(~isnan(cp))./length(cpTime),sum(~isnan(bbp))./length(bbpTime)];
fprintf('%2.1f %% of a spectra and %2.1f %% of c spectra survived\n\n',rateGood(1)*100,rateGood(2)*100)



% get times from IFCB so that sizes are the same



%% Median filter 
winSizeIn = 30; % Number of minutes on each side of time bin center
winSize   = winSizeIn/(24*60);
% temp      = datevec(min([cpTime;bbpTime]));
% startTime = datenum([temp(1:4),0,0]);
% temp      = datevec(max([cpTime;bbpTime]));
% endTime   = datenum([temp(1:3),temp(4)+1,0,0]);
% timeHour  = [startTime + 0.5/24:(1/24):endTime + 0.5/24]; % Bin center on :30

% use ifcb times
temp           = datevec(min(ifcbTime));
startTime      = datenum([temp(1:4),0,0]);
temp           = datevec(max(ifcbTime));
endTime        = datenum([temp(1:3),temp(4)+1,0,0]);
timeHour       = [startTime:(1/24):endTime];


binData      = NaN(length(timeHour)-1,3); % to match the way it was done for the ifcb
binData(:,1) = timeHour(1:end-1)';
for iHour = 1:length(timeHour)-1 % to match the way it was done for the ifcb
    indCp  = (cpTime >= timeHour(iHour)-winSize) & (cpTime <= timeHour(iHour)+winSize);
    indBbp = (bbpTime >= timeHour(iHour)-winSize) & (bbpTime <= timeHour(iHour)+winSize);
    binData(iHour,2) = nanmedian(cp(indCp));
    binData(iHour,3) = nanmedian(bbp(indBbp));
end

% Updated Figure
figsiz = [0 0 10 6].*1.6;
figure('Units','inches','Position',figsiz,...
    'PaperSize',figsiz(3:4),'PaperPosition',figsiz);
%tlo = tiledlayout(2,1,'TileSpacing','none','Padding','compact');
%nexttile;
yyaxis left
plot(binData(:,1),binData(:,2),'.-')
%nexttile;
yyaxis right
plot(binData(:,1),binData(:,3),'.-')

cd(procdata)

cstarcp = binData(:,2);
ecobbp = (binData(:,3)).*2*pi*1.077; 
gmt = binData(:,1);

save cp_bbp_hourly.mat cstarcp ecobbp gmt