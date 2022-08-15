function [Qc] = get_Qc(diam_um,realn,nprime,watern,lambda_um)

M = complex(realn,nprime);
refwater=1.3308; % Water refractive index
lambda=lambda*10^-9; % cstar wavelength (in meters)

diameter = diam_um.*10^-6; % to meters
%% particle radius
r = (diameter./2); % m

%% calculate scalar size parameter
X=2*pi*r*refwater/lambda; %dimensionless size parameter, accounting for water index of refraction

% calculate attenuation efficiencies (Qc) - only care about Qc below for
% now:
for i = 1:length(X)
    [S1(:,i),S2(:,i),Qb(i),Qc(i),Qbb(i)] = fastmie(X(i),M,90) ;
end
