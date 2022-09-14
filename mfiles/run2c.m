function run2c(ifcbdata,slope,No,refwater,nimag,lambda,refdiam,cstarcp,lutfile,nfile)

%% RUN 1 
load(ifcbdata,'bins','binwidth','ifcbPsdnBinHr','ifcbTimeHr')
data = ifcbPsdnBinHr; 
gmt = ifcbTimeHr; clear ifcbPsdnBinHr ifcbTimeHr

% build full PSD:
% first, clear bins 1:18 for clarity (we don't trust IFCB data at those
% bins)
data(:,[1:refdiam-1])=NaN;

% get mean diameter of extrapolated portion of psd. https://www.oceanopticsbook.info/view/optical-constituents-of-the-ocean/level-3/creating-particle-size-distributions-data
dias = ((1-slope)*(bins(:,2).^(2-slope) - bins(:,1).^(2-slope))) ./ ((2-slope)*(bins(:,2).^(1-slope) - bins(:,1).^(1-slope)));
dias = dias';

% extrapolate differently for 0.2-2um and 2-4um size fractions.
% assume 0.2:2um size fraction does not change over time, so assume a fixed
% slope and No
% assume 2-4um size fraction oscillates over the diel proportionally to the
% 4um size class.
N02_2 = No.*((dias(1:14)./dias(refdiam)).^-slope);
N2_6 = data(:,refdiam).*((dias([15:refdiam-1])./dias(refdiam)).^-slope);

% add those extrapolated portions to the data to make a full psd
fullpsdn = data;
fullpsdn(:,[1:14])=repmat(N02_2,length(fullpsdn(:,1)),1);
fullpsdn(:,[15:refdiam-1])=N2_6;

% add same gaps as measured data jsut in case this biases things later on
ind = isnan(fullpsdn(:,refdiam));
fullpsdn(ind,[1:refdiam-1])=NaN;
fullpsd = fullpsdn.*binwidth'; clear data fullpsdn

% now, ready to calculate cpp for the fullpsdn generated above, for a range
% of indices of refraction. This will be the LUT to be compared against
% measured cp to determine the n that gives the closest cp result to the
% measured cp.
[cp Qc Qcstar nrange dias] = calc_cp(fullpsd,dias,nimag,lambda,refwater); % saves also beamc, gmt, cyc/acyc indices, ifcb...will have to improve this in the future

save(lutfile,'cp','Qc','Qcstar','nrange','dias','slope','No','gmt');
disp('LUT file saved')

% estimate n
[n_ideal ~] = Get_Estimated_n(cstarcp,cp,nrange);

gmt = gmt(:);
save(nfile,'n_ideal','gmt','cstarcp','slope','No')


end