%% 

mdir=cd;
cd ..
cd data
datadir = cd;
cd ..\proc
procdata = cd;
cd(mdir)

load([procdata '\cyc_acyc_ind.mat'])
load([datadir '\all_features_small.mat']);
ifcb = data; clear data
ifcbTime    = datenum(ifcb.Date);
ifcbVolAn   = ifcb.Volume_analyzed.*10^-3; % cm3 to L
ifcbBiovol  = ifcb.Biovolume.*10^-9; % um3 to uL 
ifcbDiam    = ifcb.ESD; %um, from BIOVOLUME
%ifcbDiam    = (ifcb.CSA./pi).^(1/2); %um, from Cross-sectional area
ifcbVolConc = ifcbBiovol./ifcbVolAn; % ul/L --> Accounts for Volume Analyzed


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


cyc_ind = find(ifcbTime>=gmt(cyc(1))&ifcbTime<=gmt(cyc(length(cyc))));
acyc_ind = find(ifcbTime>=gmt(acyc(1))&ifcbTime<=gmt(acyc(length(acyc))));


% Get bin sizes
x         = 0.2023; % lower limit 
rho       = 200^(1/32);
bins(:,1) = x*rho.^([0:37]); % lower limit
bins(:,2) = x*rho.^([1:38]); % upper limit
binwidth  = bins(:,2) - bins(:,1);

% all:
ifcbPsdFull     = NaN(1,38);

[it,uu]        = unique(ifcbTime);
volanalyzedFull = sum(ifcbVolAn(uu));
for iBin = 1:38 % good data start at 4.3um
    [ind,indd] = find((ifcbDiam >= bins(iBin,1)) & (ifcbDiam < bins(iBin,2)));
    ifcbPsdFull(1,iBin)           = length(ind)./volanalyzedFull; % #/L, account for vol analyzed
    ifcbPsdCountFull(1,iBin)           = length(ind) % #
end
ifcbPsdnFull = ifcbPsdFull./binwidth';


% cyclne
ifcbPsdFull_cyc     = NaN(1,38);
ifcbTime_cyc = ifcbTime(cyc_ind);
ifcbVolAn_cyc = ifcbVolAn(cyc_ind);
ifcbDiam_cyc = ifcbDiam(cyc_ind);

[it,uu]        = unique(ifcbTime_cyc);
volanalyzedFull = sum(ifcbVolAn_cyc(uu));
for iBin = 1:38 % good data start at 4.3um
    [ind,indd] = find((ifcbDiam_cyc >= bins(iBin,1)) & (ifcbDiam_cyc < bins(iBin,2)));
    ifcbPsdFull_cyc(1,iBin)           = length(ind)./volanalyzedFull; % #/L, account for vol analyzed
    ifcbCountFull_cyc(1,iBin)                   = length(ind); % # per bin ,get counts for stats
end
ifcbPsdnFull_cyc = ifcbPsdFull_cyc./binwidth';

% anticyclone
ifcbPsdFull_acyc     = NaN(1,38);
ifcbTime_acyc = ifcbTime(acyc_ind);
ifcbVolAn_acyc = ifcbVolAn(acyc_ind);
ifcbDiam_acyc = ifcbDiam(acyc_ind);

[it,uu]        = unique(ifcbTime_acyc);
volanalyzedFull = sum(ifcbVolAn_acyc(uu));
for iBin = 1:38 % good data start at 4.3um
    [ind,indd] = find((ifcbDiam_acyc >= bins(iBin,1)) & (ifcbDiam_acyc < bins(iBin,2)));
    ifcbPsdFull_acyc(1,iBin)           = length(ind)./volanalyzedFull; % #/L, account for vol analyzed
    ifcbCountFull_acyc(1,iBin)                   = length(ind); % # per bin ,get counts for stats
end
ifcbPsdnFull_acyc = ifcbPsdFull_acyc./binwidth';





% srise vs sunset PSD
t=ifcbTime-10/24;
tt = datevec(t);
hrise = find(tt(:,4)==7);
hset = find(tt(:,4)==20);

% srise
ifcbPsdFull_srise    = NaN(1,38);
ifcbTime_srise = ifcbTime(hrise);
ifcbVolAn_srise = ifcbVolAn(hrise);
ifcbDiam_srise = ifcbDiam(hrise);


[it,uu]        = unique(ifcbTime_srise);
volanalyzedFull = sum(ifcbVolAn_srise(uu));
for iBin = 1:38 % good data start at 4.3um
    [ind,indd] = find((ifcbDiam_srise >= bins(iBin,1)) & (ifcbDiam_srise < bins(iBin,2)));
    ifcbPsdFull_srise(1,iBin)           = length(ind)./volanalyzedFull; % #/L, account for vol analyzed
    ifcbCountFull_srise(1,iBin)                   = length(ind); % # per bin ,get counts for stats
end
ifcbPsdnFull_srise = ifcbPsdFull_srise./binwidth';

% sunset
ifcbPsdFull_sset    = NaN(1,38);
ifcbTime_sset = ifcbTime(hset);
ifcbVolAn_sset = ifcbVolAn(hset);
ifcbDiam_sset = ifcbDiam(hset);

[it,uu]        = unique(ifcbTime_sset);
volanalyzedFull = sum(ifcbVolAn_sset(uu));
for iBin = 1:38 % good data start at 4.3um
    [ind,indd] = find((ifcbDiam_sset >= bins(iBin,1)) & (ifcbDiam_sset < bins(iBin,2)));
    ifcbPsdFull_sset(1,iBin)           = length(ind)./volanalyzedFull; % #/L, account for vol analyzed
    ifcbCountFull_sset(1,iBin)                   = length(ind); % # per bin ,get counts for stats
end
ifcbPsdnFull_sset = ifcbPsdFull_sset./binwidth';



cd(mdir)

clearvars -except ifcb* datadir mdir procdata figsdir cyc acyc