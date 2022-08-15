function run_1(ifcbdata,slope,No,refwater,nimag,lambda,refdias,lutfile,nfile)

%% RUN 1 
load(ifcbdata,'bins','binwidth','ifcbPsdnBinHr','dias','ifcbTimeHr')
data = ifcbPsdnBinHr; 
gmt = ifcbTimeHr; clear ifcbPsdnBinHr ifcbTimeHr

% build full PSD:
% first, clear bins 1:18 for clarity (we don't trust IFCB data at those
% bins)
data(:,[1:18])=NaN;

% create extrapolated portion of the PSD absed on above slope and No
N = No.*((dias(1:18)./dias(refdias)).^-slope);

% assume PSD at bins 1:18 do not change over time, so copy same value for
% all rows
fullpsdn = data;
fullpsdn(:,[1:18])=repmat(N,length(fullpsdn(:,1)),1);
% add same gaps as measured data jsut in case this biases things later on
ind = isnan(fullpsdn(:,19));
fullpsdn(ind,[1:18])=NaN;
fullpsd = fullpsdn.*binwidth'; clear data fullpsdn


% now, ready to calculate cpp for the fullpsdn generated above, for a range
% of indices of refraction. This will be the LUT to be compared against
% measured cp to determine the n that gives the closest cp result to the
% measured cp.
[cp Qc Qcstar nrange dias] = calc_cp(fullpsd,bins,slope,nimag,lambda,refwater); % saves also beamc, gmt, cyc/acyc indices, ifcb...will have to improve this in the future


save(lutfile,'cp','Qc','Qcstar','nrange','dias','slope','No','gmt')

% estimate n
load('cp_bbp_hourly.mat')
[n_ideal] = Get_Estimated_n(cstarcp,cp,nrange);

gmt = gmt(:);
save(nfile,'n_ideal','gmt','cstarcp','ecobbp','slope','No')
end