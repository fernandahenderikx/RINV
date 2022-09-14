% FASTMIE Vectorized version of BHMIE 
%   Calculates the amplitude scattering matrix elements and efficiencies
%   for extinction, total scattering, and backscattering for a given size
%   parameter and relative refractive index. This code maintains the
%   calling syntax of E. Boss' BHMIE translation.
%
%   [S1,S2,Qb,Qc,Qbb] = FASTMIE(X,M,NANG)
%   X and M are the scalar size parameter (2*pi*r/lambda) and relative
%   refractive index of the particle. NANG is the number of angles between
%   0 and 90 degrees; matrix elements are calculated at 2*NANG-1 angles
%   including 0, 90, and 180 degrees.
%
%   [S1,S2,Qb,Qc,Qbb] = FASTMIE(X,M,[],THETA)
%   Same as above, but use user-specified THETA vector (radians).
%
%   FASTMIE returns
%   S1 and S2 are 2*NANG-1 x 1 vectors of scattering matrix elements;
%   Qb, Qc, and Qbb are scalar extinction, total scattering, and
%   backscattering efficiencies.
%
% Wayne H. Slade
% University of Maine
% 10 Aug 2005

% 3 May 2006 added ability to specify angles for which to calculate
% scattering matrix elements

% 31 May 2022 Added in P11 for VSF and acceptance angle calculations


function [s1, s2, Q_scat, Q_ext, Q_back, Q_extStar] = fastmie_jga(x,m,nang,varargin)

if any([length(m) length(x)] ~= 1)
    error('m, x must be scalar')
end

if ~(length(nang) == 1 || isempty(nang))
    error('nang must be scalar or empty matrix')
end

if isempty(nang)
    theta = varargin{1}(:);
else
    theta = linspace(0,pi,nang*2 - 1); % calculated at (nang*2-1) angles for agreement with bhmie
end

nstop = ceil(2 + x + 4*x.^(1/3));
n     = (1:nstop)';

jx  = besselj(n + 0.5,x)  .*sqrt(0.5*pi/x);
jmx = besselj(n + 0.5,m*x).*sqrt(0.5*pi/(m*x));
yx  = bessely(n + 0.5,x)  .*sqrt(0.5*pi/x);
hx  = jx + j*yx;

j1x  = [sin(x)/x;       jx(1:nstop - 1)];
j1mx = [sin(m*x)/(m*x); jmx(1:nstop - 1)];
y1x  = [-cos(x)/x;      yx(1:nstop - 1)];
h1x  = j1x + j*y1x;

d_jx  = x.*j1x - n.*jx;          % derivative of spherical bessel j_n(x)
d_jmx = (m*x).*j1mx - n.*jmx;    % derivative of spherical bessel j_n(mx)
d_hx  = x.*h1x - n.*hx;

an = ((m^2).*jmx.*d_jx - jx.*d_jmx)./((m^2).*jmx.*d_hx - hx.*d_jmx);
bn = (jmx.*d_jx - jx.*d_jmx)./(jmx.*d_hx - hx.*d_jmx);

% calculate pi and tau
[pi_n,tau_n] = mie_pt(theta,nstop);   % pi and tau will be (nmax x length(theta))

p  = repmat((2*n + 1)./(n.*(n + 1)),1,length(theta)).*pi_n;
t  = repmat((2*n + 1)./(n.*(n + 1)),1,length(theta)).*tau_n;
AN = repmat(an,1,length(theta));
BN = repmat(bn,1,length(theta));
s1 = sum(AN.*p + BN.*t).'; % transpose s1,s2 for dim agreement with bhmie
s2 = sum(AN.*t + BN.*p).';

Q_scat = (2/x.^2) * sum((2*n + 1).*(abs(an).^2 + abs(bn).^2)); % bh p103
Q_ext  = (2/x.^2) * sum((2*n + 1).*real(an + bn));
Q_back = (1/x.^2) * abs(sum((2*n + 1).*((-1).^n).*(an - bn))).^2;

% Calculate P11 of Mueller Matrix
P11        = 0.5*(abs(s1).^2 + abs(s2).^2)';
ang        = theta;
angDeg     = ang/pi*180;
ntgrand    = P11.*sin(ang);
P11IntFull = trapz(angDeg,ntgrand,2);
indAcc     = find(angDeg >= 1.2);
P11IntAcc  = trapz(angDeg(indAcc),ntgrand(indAcc),2);
Q_extStar  = Q_ext.*(P11IntAcc/P11IntFull);

return

function [pi_n,tau_n] = mie_pt(theta,nmax)
% Bohren and Huffman (1983), p. 94 - 95
% (adapted the code from Christian Matzler)
%
% theta (radians) should be (1 x M), where M is the number of angles to consider
% nmax is scalar
% pi and tau will be (nmax x M)

theta = theta(:)';   % ensure that cos_theta is a row vector
u     = cos(theta);

pi_n       = zeros(nmax,length(u));
tau_n      = zeros(nmax,length(u));
pi_n(1,:)  = 1;
tau_n(1,:) = u;
pi_n(2,:)  = 3*u;
tau_n(2,:) = 3*cos(2*acos(u));

for n = 3:nmax
    p1         = (2*n - 1)./(n - 1).*pi_n(n - 1,:).*u;
    p2         = n./(n - 1).*pi_n(n - 2,:);
    pi_n(n,:)  = p1 - p2;
    t1         = n.*u.*pi_n(n,:);
    t2         = (n + 1).*pi_n(n - 1,:);
    tau_n(n,:) = t1 - t2;
end
return



