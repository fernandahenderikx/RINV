     
% READPSCAT Read uncalibrated scattered power data for LISST-B instrument
%
% This function reads corrected (but uncalibrated) scattered power data for
% each detector ring, along with beam attenuation.
%
% [ANGR, PSCAT, BEAMC, SNR] = READPSCAT(DATAFILE, ZSCFILE, DCALFILE)
%   -or- [...] = READPSCAT(..., BADSCANS, BADRINGS, REJECTSNR)
%   -or- [ANGR,PSCAT,BEAMC,SNR, ANG1,ANG2,RINGS] = READPSCAT(...) 
%
% Inputs are the data file to be read, zscat file (see notes), and the 
% ringarea file.
%
% Optionally, Vectors of scan numbers and ring detectors to reject can be
% supplied, for example:
%   [...] = READPSCAT(DATAFILE,ZSCFILE,DCALFILE, [50:100], [1:6], 50)
% 
% The SNR threshold for rejecting a ring can also be specified, otherwise
% all rings (not in BADRINGS) are kept. REJECTSNR is not in dB; SNR is 
% calculated for raw data as mean(data)/std(zscat).
%
% The ZSCFILE specified can be a Sequoia-style single average of scans, or 
% can be a log file taken on clean water (not corrected for any zscat). It's
% easy: use a zscat file loaded with zeros and log data as you normally 
% would with the LISST chamber full of particle free DIW! If a 
% log file is specified, the std of that data is used in calculating the
% snr, otherwise a set of default std(zscat) are used. 
%
% The function returns:
%   ANGR (1 x nring), a vector of mean ring angles in radians
%   PSCAT (nscan x nring) a matrix of corrected scatter values for the
%         nring ring detectors not rejected over nscan scans
%   BEAMC (nscan x 1) contains measured values of beam attenuation
%   SNR (1 x nring) a vector of signal to noise ratios for the rings that
%       were not discarded
%   ANG1,ANG2 each (1 x nring) giving the detector edge angles in radians
%   RINGS (1 x nring) giving the ring numbers for the rings that are kept
% 
%
% Wayne H. Slade
% University of Maine, School of Marine Sciences
% Contact wayne.slade@gmail.com

% 10 mar 06 not removing scans with negative data
% 15 mar 06 modified code such that angles returned are IN WATER
% 24 Sept 2021 modified code to output tau and raw data (JGA)

function [angr,pscat,beamc,snr,ang1,ang2,rings,cscat,tau,rawdata] = ...
    readpscat(datafile, zscfile, dcalfile, varargin)

	%disp(['Reading DATAFILE: ' datafile])
	%disp(['         ZSCFILE: ' zscfile])
	%disp(['        DCALFILE: ' dcalfile])

	% Handle variable function arguments
	error(nargchk(3,6, nargin, 'struct'))
	if nargin == 4
		badscans  = varargin{1};
		badrings  = [];
		rejectsnr = [];
	elseif nargin == 5
		badscans  = varargin{1};
		badrings  = varargin{2};
		rejectsnr = [];
	elseif nargin == 6
		badscans  = varargin{1};
		badrings  = varargin{2};
		rejectsnr = varargin{3};
	else
		badscans  = [];
		badrings  = [];
		rejectsnr = [];
	end
	
	% Logarithmically spaced edges of ring detectors (radian)
	%theta_i = logspace(0,log10(200),33)*(0.097/60);
	%theta_i = logspace(0,log10(200),33)*(0.1/60);
	theta0air = 0.1/180*pi;
	theta_i   = (200.^((0:32)/32)).*theta0air;
	theta_j   = theta_i(2:33);
	theta_i   = theta_i(1:32);
	
	theta_i = asin(sin(theta_i)./1.3308); 
	theta_j = asin(sin(theta_j)./1.3308); 
	
	% Fraction of circle covered by detector
	phi = 1/6;
	
    % Instrument pathlength [m]
	pathl = 0.05;
	rings = 1:32;
	
	% Load data files-- 
	% uncalibrated scatter (counts) for each ring (won't look like vsf)
	data = load(datafile);

	% "ringarea" file; actually just correction factors for detector
	% responsivity
	dcal = load(dcalfile);
	% uncalibrated scatter for each ring for clean particle-free water
	% note that the zscat should be logged as a series of scans in order to
	% determine the variability of the measurement system
	zscat = load(zscfile);
	
	% Ensure that data is [nscan x 40]
	if size(data,2)~=40
		error('DATA not [nscan x 40]')
	end
	nscan   = size(data,1);
    rawdata = data;
	%disp(['DATAFILE contains ' num2str(nscan) ' scans...'])
	
	% Ensure that dcal is 32 element row vector
	if size(dcal,2) ~= 32
		error('DCAL not [nscan x 32]')
	end
	
	% Ensure that zscat has 40 columns
	if ~(all(size(zscat) == [40 1]) || size(zscat,2) == 40)
		error('ZSCAT not [nscan x 40] or [40 x 1]')
	end
	
	% Two possibilities for zscat: (1) Sequoia style single [1x40] average
	% of several scans, or (2) clean water data file [m x 40]
	% If (1) then snr is based on a default* std(zsc), otherwise use the
	% real data to calculate std(zsc)
    
    % *Update from James: I hate this, cause you have to manually load in a
    % datafile where you run clean water instead of a zscat. Luckily, I did
    % this during the cruise... Need to add to protocol. 
	if all(size(zscat) == [40 1])
		% We're assuming that variability in the 1KDalton comparison file is 
		%   representative of instrument noise in clean water.
		%disp('!! USING DEFAULT STD(ZSC) for S/N 1285 !!')
		stdzsc = [2.6508,0.6294,0.4088,0.2083,0.2420,0.7042,0.5771,0.2177,...
            0.4071,0.3046,0.1902,0.2818,0.2206,0.1921,0.1746,0.1671,0.2230,...
            0.1496,0.1723,0.1283,0.1689,0.1428,0.1718,0.1207,0.1365,0.0985,...
            0.1009,0.0989,0.0995,0.0865,0.0985,0.1216];
		%stdzsc = repmat(0.500, 1,32);
		
		% And make sure that zscat is a row vector instead of column
		zscat = zscat(:)';
	else
		%disp('Using ZSCAT data to determine STD(ZSC)...')
		stdzsc = std(zscat(:,1:32));
		
		% And make zscat a row vector, mean of all scans
		zscat = mean(zscat,1);
	end				  
	
	% Sometimes rings (large angles usually) have sero std, set their std
	% to nan
	stdzsc(stdzsc==0) = nan;
	
	% A reasonable reject criteria is SNR < ~40 which is approximately 32dB
	% (see run_snrtest.m). 
	snr = mean(data(:,1:32)) ./ stdzsc;
	% snr = mean(data(:,1:32)) ./ std(data(:,1:32));
	% snr = mean(data(:,1:32)) ./ max([std(zscat(:,1:32)); std(data(:,1:32))]);

	% For more discussion of the following data processing, see Sequoia 
	% App Notes 6 and L012 ...

	% ratio of water-transmitted to reference laser counts in clean water
	rr = zscat(33)/zscat(36);
	% transmission and beam attenuation
	tau = data(:,33)./rr./data(:,36);
	beamc = -log(tau)/pathl;

	% "calibrated scatter" implies that we've corrected the raw scatter for
	% the detector area (determined based on geometry), zscat, dcal
	% (responsivity correction, aka ringarea file), and attenuation 
	% within the sample volume:
	% (1) first correct for attenuation (using tau) and subtract zscat...
	%     code in MatlabBundle.zip/nlia_func.m also has a z(36)/data(1,36)
	%     factor to correct for drift in laser power; we make the 
	%     correction later similar to vsfcode_new.m
	cscat = data(:,1:32)./repmat(tau,1,32) - repmat(zscat(1:32),nscan,1);
	% (2) second correct for detector responsivities
	cscat = cscat .* repmat(dcal,nscan,1);
	% (3) correct for laser power (counts) in water (relative to zscat)
	cscat = cscat ./ repmat(data(:,36)/mean(zscat(36)),1,32);
	% (4) apply analytic formulation for ring area
	pscat = cscat ./ repmat(pi*phi*pathl*(theta_j.^2-theta_i.^2),nscan,1);
	
	% mean detector angles (radian)
    angr = mean([theta_i; theta_j]);
	% detector edge angles
	ang1 = theta_i;
	ang2 = theta_j;
	
	% Remove data not meating quality criteria or corresponding to BADRINGS
	% or BADSCANS --
	
	% (1) discard data scans that are specified as BADSCANS
	if ~isempty(badscans)
		pscat(badscans,:) = [];
		beamc(badscans,:) = [];
		disp(['** Discard BADSCANS: ' num2str(badscans)])
	end
	
    % (2) discard rings that are specified as BADRINGS
	if ~isempty(badrings)
		pscat(:,badrings) = [];
		angr(:,badrings) = [];
		ang1(:,badrings) = [];
		ang2(:,badrings) = [];
		rings(:,badrings) = [];
		snr(:,badrings) = [];   
		disp(['** Discard BADRINGS: ' num2str(badscans)])
	end
	
	% (3) reject detectors based on low SNR criteria
	if ~isempty(rejectsnr)
		lowsnr = snr < rejectsnr;
		pscat(:,lowsnr) = [];
		angr(lowsnr) = [];
		ang1(lowsnr) = [];
		ang2(lowsnr) = [];
		rings(lowsnr) = [];
		snr(lowsnr) = [];
		disp(['** Discard due to low snr rings: ' num2str(find(lowsnr))])
	end
	
	% % (4) discard any remaining scans with negative in cscat or cmeas
	% badscan = any(cscat'<=0)' | cmeas<=0;
	% cscat(badscan,:) = [];
	% beamc(badscan) = [];
	% disp(['** Discard due to negative in scan: ' num2str(length(find(badscan))) ' rows'])
	
	