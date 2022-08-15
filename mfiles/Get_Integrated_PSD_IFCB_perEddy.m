%% 

mdir=cd;
cd ..
cd data
datadir = cd;
cd ..\proc
procdata = cd;
cd(mdir)

load([procdata '\cyc_acyc_ind.mat'])
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

cyc_ind = find(ifcbTime>=gmt(cyc(1))&ifcbTime<=gmt(cyc(length(cyc))));
acyc_ind = find(ifcbTime>=gmt(acyc(1))&ifcbTime<=gmt(acyc(length(acyc))));


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

% all:
ifcbPsdFull     = NaN(1,38);

[it,uu]        = unique(ifcbTime);
volanalyzedFull = sum(ifcbVolAn(uu));
for iBin = 1:38 % good data start at 4.3um
    [ind,indd] = find((ifcbDiam >= bins(iBin,1)) & (ifcbDiam < bins(iBin,2)));
    ifcbPsdFull(1,iBin)           = length(ind)./volanalyzedFull; % #/L, account for vol analyzed
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
end
ifcbPsdnFull_sset = ifcbPsdFull_sset./binwidth';



cd(mdir)

clearvars -except ifcbPsdnFull* datadir mdir procdata figsdir cyc acyc dias