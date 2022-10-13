
clear all; close all;
mdir=cd;
cd ..\..
addpath(genpath('RINV'))
cd(mdir)
cd ..
cd data
datadir = cd;
cd ..\proc
procdata = cd;
cd ..\figures
figsdir = cd;
cd(mdir)

cd(procdata)
load('ifcb_cruiseIntegrated.mat','bins','ifcbPsdnFull','ifcbPsdFull')

slope = 4.57;
No = 5.3113e+04;
refdiam = 21;

% get mean diameter of extrapolated portion of psd. https://www.oceanopticsbook.info/view/optical-constituents-of-the-ocean/level-3/creating-particle-size-distributions-data
dias = ((1-slope)*(bins(:,2).^(2-slope) - bins(:,1).^(2-slope))) ./ ((2-slope)*(bins(:,2).^(1-slope) - bins(:,1).^(1-slope)));
dias = dias';
refdiam_value = dias(refdiam);

N02_6 = No.*((dias(1:refdiam-1)./dias(refdiam)).^-slope);

% add those extrapolated portions to the data to make a full psd
fullpsdn = ifcbPsdnFull;
fullpsdn(:,[1:refdiam-1])=repmat(N02_6,length(fullpsdn(:,1)),1);

save refdiam.mat dias fullpsdn refdiam refdiam_value slope No