%% Create cp LUT for KM1709
%% calculate cp from IFCB (measured+extrapolated) PSD and for a range of n

function [cp Qc_all Qcstar_all nall dias] = calc_cp(fullpsd,dias,nimag,lambda,refwater)

% get lambda and diameter in meters
lambda=lambda*10^-9; % CSTAR wavelength (in meters)
nint = 0.0001; % n spacing
diameter= dias*10^-6; % meters;;


nall = [1.001:nint:1.2]; % range o n to run analysis for
%nall = [1.02:nint:1.04]; % limited range to speed it up when needed for testing purposes

psddata = fullpsd;


clear cp Qc Qcstar

for nn = 1:length(nall)

        %% calculate complex relative refractive index of the particle
    M = complex(nall(nn),nimag);

    %% particle radius
    r = (diameter./2); % m

    %% calculate scalar size parameter
    X=2*pi*r*refwater/lambda; %dimensionless size parameter, accounting for water index of refraction

    % calculate attenuation efficiencies (Qc) - only care about Qc below for
    % now:
        for i = 1:length(X)
           [~,~,~,Qc(i),~,Qcstar(i),~] = fastmie_mod(X(i),M,3601); %modifiend from fastmie.m to calculate Qc* and Q11
        end
    clear S1 S2 Qb Qbb

    %% cross section * Qc
    cs = pi.*(r.^2).*Qcstar; % m2

    %% get psd in #/m^3
    N = psddata.*1000;

    %% cp per bin
    cp_perbin= N.*cs; % m-1

    %% get total cp
        for j = 1:length(psddata(:,1))
            cp(j,nn)=nansum(cp_perbin(j,:)); %
        end
        disp(nall(nn))
        

        Qc_all(nn,:)=Qc;
        Qcstar_all(nn,:)=Qcstar;

end
        cp(cp==0)=NaN;
end

