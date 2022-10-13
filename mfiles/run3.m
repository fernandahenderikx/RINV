function run3(ifcbdata,slope,No,refwater,nimag,lambda,refdiam,cstarcp,lutfile,nfile)

%% RUN 3 (truncate PSD to 10um, 20 um, 50 um) 
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
N02_6 = No.*((dias(1:refdiam-1)./dias(refdiam)).^-slope);

% add those extrapolated portions to the data to make a full psd
fullpsdn = data;
fullpsdn(:,[1:refdiam-1])=repmat(N02_6,length(fullpsdn(:,1)),1);

% add same gaps as measured data jsut in case this biases things later on
ind = isnan(fullpsdn(:,refdiam));
fullpsdn(ind,[1:refdiam-1])=NaN;
fullpsd = fullpsdn.*binwidth'; clear data fullpsdn




% truncate data at 10 um, 20um, 50 um to see how the estimated n varies
% when large particles are ignored (or not)
for fff = 1:3 % 1:4
    if fff == 1
        fullpsd_use = fullpsd(:,[1:24]);
        dias_use = dias([1:24]);
        nfile_use = [nfile(1:end-4) '_max10um.mat'];
        lutfile_use = [lutfile(1:end-4) '_max10um.mat'];
    end
    if fff == 2
        fullpsd_use = fullpsd(:,[1:28]);
        dias_use = dias([1:28]);
        nfile_use = [nfile(1:end-4) '_max20um.mat'];
        lutfile_use = [lutfile(1:end-4) '_max20um.mat'];
    end
    if fff == 3
        fullpsd_use = fullpsd(:,[1:34]);
        dias_use = dias([1:34]);
        nfile_use = [nfile(1:end-4) '_max50um.mat'];
        lutfile_use = [lutfile(1:end-4) '_max50um.mat'];
    end
    if fff == 4
        fullpsd_use = fullpsd(:,[1:20]);
        dias_use = dias([1:20]);
        nfile_use = [nfile(1:end-4) '_max6um.mat'];
        lutfile_use = [lutfile(1:end-4) '_max6um.mat'];
    end

% now, ready to calculate cpp for the fullpsdn generated above, for a range
% of indices of refraction. This will be the LUT to be compared against
% measured cp to determine the n that gives the closest cp result to the
% measured cp.
[cp Qc Qcstar nrange diam] = calc_cp(fullpsd_use,dias_use,nimag,lambda,refwater); % saves also beamc, gmt, cyc/acyc indices, ifcb...will have to improve this in the future
    
    save(lutfile_use,'cp','Qc','Qcstar','nrange','diam','slope','No','gmt')
    
    % estimate n
    [n_ideal ~] = Get_Estimated_n(cstarcp,cp,nrange);
    
    gmt = gmt(:);
    save(nfile_use,'n_ideal','gmt','cstarcp','slope','No')
    end
end