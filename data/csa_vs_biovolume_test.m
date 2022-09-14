%% check diameter estimated from CSA vs diameter estimated from biovolume to see that as diameter gets bigger the two estimates deviate the most

load all_features_small.mat

d_csa = (data.CSA./pi).^(1/2).*2;
d_vol = (data.Biovolume./(pi*3/4)).^(1/3).*2;

figure
subplot 121
plot(d_vol,d_csa,'.')
xlabel('diameter estimated from cross-sectional area')
ylabel('diameter estimated from volume')

subplot 122
hist(d_vol-d_csa,1000)
xlim([-5 0.5])
title('diameter from volume - diameter from csa')