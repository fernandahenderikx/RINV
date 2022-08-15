function nsw = get_water_n(lambda,Tc,S)
% Test Variables
% lambda = 660; % in nm
% Tc     = 10; % in degrees
% S      = 33; % in PSU
% 
% refractive index of air is from Ciddor (1996,Applied Optics)
n_air = 1.0 + (5792105.0./(238.0185 - 1./(lambda/1e3).^2) + ...
        167917.0./(57.362 - 1./(lambda/1e3).^2))/1e8;

% refractive index of seawater is from Quan and Fry (1994, Applied Optics)
n0 = 1.31405; 
n1 = 1.779e-4; 
n2 = -1.05e-6; 
n3 = 1.6e-8; 
n4 = -2.02e-6;
n5 = 15.868; 
n6 = 0.01155;  
n7 = -0.00423;  
n8 = -4382; 
n9 = 1.1455e6;

nsw = n0 + (n1 + n2*Tc + n3*Tc.^2).*S + n4*Tc.^2 + (n5 + n6*S + n7*Tc)./lambda + ...
      n8./lambda.^2 + n9./lambda.^3; % pure seawater
nsw = nsw.*n_air;
