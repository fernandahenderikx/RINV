%% MAIN script for Index of Refraction paper

% main directory should contain the following folders: mfiles, data, proc,
% figures.

clear all; close all;
mdir=cd;
cd ..
cd data
datadir = cd;
cd ..\proc
procdata = cd;
cd ..\figures
figsdir = cd;
cd(mdir)


%% load and process min avgd beamc and bbp, cleanup, make hourly
cd(mdir)
CleanupIOPs % get binData with hourly time, cp, bbp

%% process IFCB data (integrated hourly averages for analysis of n)
cd(mdir)
Get_HourlyIntegrated_PSD_IFCB % get 415x38 PSDs, Volconcentration, Median diameters, etc.

%% process IFCB data (Integrated cruise averages to get PSD)
cd(mdir)
Get_CruiseIntegrated_PSD_IFCB % generates ifcbPSDfitReport.txt report; also report generated at the end of script which shows slope and No for uncertainty analysis of the extrapolated portion of the PSD

% Slope: 4.7394, SE: 0.1182, STD: 0.4094
% N0:10^4.9046, SE of exponent: 0.0552, STD of exponent: 0.1912
% 
% N0(+1sigma): 1.246765e+05
% N0: 8.028211e+04
% N0(-1sigma): 5.169553e+04


%%  Estimate n based on different assumptions for the extrapolated portion of the IFCB
cd(procdata)
tic % this takes ~5hours after changing nangle=3601 to calculate Qc* seconds

% fixed parameters
ifcbdata = 'ifcb_IntegratedHourly.mat';
%refwater = get_water_n(660,26,35) % == 1.337
refwater=1.337; % Water refractive index
lambda = 660; % nm; will convert to meter inside code below
refdias = 19; % reference diameter. This was used to determine No below, so do not change that.
% % choice of n':
% lambdaair=lambda/1.337; % in water
% nimag= 0.010658*exp(-0.007186*lambda) %  == 0.0003; Stramski et al 2004
nimag = 0.0003;

% define parameters of extrapolated portion of PSD:
% this is the main result to be discussed: average psd and average No
slope = 4.7;
No = 10^4.9046;
lutfile = [procdata '\LUTcp_slope47_NoAvg.mat']; % name to save as
nfile = [procdata '\n_slope47_NoAvg.mat']; % name to save as
run1(ifcbdata,slope,No,refwater,nimag,lambda,refdias,lutfile,nfile);

% allow 2-4um to oscillate over diel
lutfile = [procdata '\LUTcp_slope47_NoAvg_diel2to4.mat']; % name to save as
nfile = [procdata '\n_slope47_NoAvg_diel2to4.mat']; % name to save as
run2(ifcbdata,slope,No,refwater,nimag,lambda,refdias,lutfile,nfile)

% vary No
slope = 4.7;
No = 1.246765e+05;
lutfile = [procdata '\LUTcp_slope47_NoPlusOneSigma.mat']; % name to save as
nfile = [procdata '\n_slope47_NoPlusOneSigma.mat']; % name to save as
run1(ifcbdata,slope,No,refwater,nimag,lambda,refdias,lutfile,nfile);

slope = 4.7;
No =5.169553e+04;
lutfile = [procdata '\LUTcp_slope47_NoMinusOneSigma.mat']; % name to save as
nfile = [procdata '\n_slope47_NoMinusOneSigma.mat']; % name to save as
run1(ifcbdata,slope,No,refwater,nimag,lambda,refdias,lutfile,nfile);

% Vary slopes
slope = 4.3;
No = 10^4.9046;
lutfile = [procdata '\LUTcp_slope43_NoAvg.mat']; % name to save as
nfile = [procdata '\n_slope43_NoAvg.mat']; % name to save as
run1(ifcbdata,slope,No,refwater,nimag,lambda,refdias,lutfile,nfile);

slope = 5.1;
No = 10^4.9046;
lutfile = [procdata '\LUTcp_slope51_NoAvg.mat']; % name to save as
nfile = [procdata '\n_slope51_NoAvg.mat']; % name to save as
run1(ifcbdata,slope,No,refwater,nimag,lambda,refdias,lutfile,nfile);
toc


%% FIGURES

close all; 

%% Fig. 1
cd(procdata)
clearvars -except datadir mdir procdata figsdir cyc acyc

load('LUTcp_slope47_NoAvg.mat')

ind1 = find(round(nrange,4)==1.0500);
ind2 = find(round(nrange,4)==1.0900);

figure
loglog(dias,Qc(ind1,:),'ko-','MarkerFaceColor','y'); hold on; 
loglog(dias,Qcstar(ind1,:),'k^-','MarkerFaceColor','y'); hold on; 
loglog(dias,Qc(ind2,:),'ko-','MarkerFaceColor','c'); hold on; 
loglog(dias,Qcstar(ind2,:),'k^-','MarkerFaceColor','c'); hold on; 
xlim([10^-1 10^2])
ylabel('{\itQ_c or Q_c^*}   [unitless]')
xlabel('Diameter [\mum]'); 
legend('{\itQ_c; n} = 1.05','{\itQ_c^*; n} = 1.05','{\itQ_c; n} = 1.09','{\itQ_c^*; n} = 1.09');
cd(figsdir)
print(gcf,'Fig1_Qc.png','-dpng','-r300');
savefig(gcf,'Fig1_Qc.fig','compact')


% figure
% imagesc(log10(dias),nrange',Qc)
% colormap(cmocean('deep'))
% h=colorbar; title(h,'Q_c [unitless]')
% xlabel('Diameter [\mum]'); ylabel('real index of refraction [n]')
% set(gca,'XTick',[-0.5 0 0.5 1 1.5 2],'XTickLabel',{'0.3','1','3','10','30','100'},'FontSize',14);
% cd(figsdir)
% print(gcf,'Fig1_Qc.png','-dpng','-r300');
% savefig(gcf,'Fig1_Qc.fig','compact')

%% Fig. 2: flowchart, see pptx

%% Fig. 3
clearvars -except datadir mdir procdata figsdir cyc acyc
cd(procdata)
Get_Integrated_PSD_IFCB_perEddy %this is similar to what was done in Get_CruiseIntegrated_PSD_IFCB, but for the different eddies.
load([procdata '\ifcb_IntegratedHourly.mat'],'ifcbPsdBinHr','*Time*','*Agg*');
load('n_slope47_NoAvg.mat')

% indices for cyclone and anticyclone eddies:
cyc = [155:245]; 
acyc = [269:361];

psdsum_ifcb = nansum(ifcbPsdBinHr(:,19:38)');psdsum_ifcb(psdsum_ifcb==0)=NaN;
%a=movmean(ifcbVolconcAggBinHr,5);

bbp_660 = ecobbp.*(660/700)^(-1.0337); % Maritoren Siegel Peterson 2002
bbpratio = bbp_660./cstarcp;

tzone = 10;
t=gmt-tzone/24;
tt=datevec(t);
y=ones(1,18).*2017; m = [6 6 6 6 7 7 7 7 7 7 7 7 7 7 7 7 7 7]; d=[27:30 1:14];
ticks = datenum(y,m,d); 
tickLabels = d;

linewidth = 1;
listOfThermalColors = cmocean('grey',7);
listOfThermalColors = listOfThermalColors(1:6,:);
r = rgb('teal');
listOfThermalColors(1,:) = r; % change first color to red


figure_vertical
ha = tight_subplot(4,1,[.025 .03],[.08 .08],[.2 .2])
axes(ha(1))
h2=plot(t,cstarcp,'.-','Color','k','LineWidth',linewidth);hold on
axis tight; 
set(gca,'XTick',ticks,'XTickLabel','','FontSize',14);
xlim([datenum(2017,6,26,12,0,0) datenum(2017,7,14,12,0,0)]);
ylim([0.032 0.057])
ylabel({'c_{p} [m^-^1]'})
box on
addlabel('a)')
vline([t(cyc(1)) t(cyc(end)) t(acyc(1)) t(acyc(end))],'k:')
text(t(cyc(1))+0.5,0.059,'cyclone')
text(t(acyc(1)),0.059,'anticyclone')
axes(ha(2))
plot(t, psdsum_ifcb,'.-k','LineWidth',linewidth)
ylabel({'IFCB';'[# L^-^1]'})
axis tight
%ylim([0.5*10^5 3*10^5])
set(gca,'XTick',ticks,'XTickLabel','','FontSize',14);
xlim([datenum(2017,6,26,12,0,0) datenum(2017,7,14,12,0,0)]);
addlabel('b)')
vline([t(cyc(1)) t(cyc(end)) t(acyc(1)) t(acyc(end))],'k:')
axes(ha(3))
plot(t, ifcbVolconcAggBinHr,'k.-','LineWidth',linewidth) ;hold on%,'LineWidth',linewidth); hold on
%plot(t, a,'-k','LineWidth',linewidth);
ylabel({'IFCB VolConc';' [\muLL^-^1]'})
ylim([0.01 0.06])
set(gca,'XTick',ticks,'XTickLabel','','FontSize',14);
xlim([datenum(2017,6,26,12,0,0) datenum(2017,7,14,12,0,0)]);
vline([t(cyc(1)) t(cyc(end)) t(acyc(1)) t(acyc(end))],'k:')
addlabel('c)')
axes(ha(4))
for k = 1 : 6
    plot(t,ifcbPsdBinHr(:,k+18),'-', 'Color', listOfThermalColors(k, :), 'LineWidth', 1.5,'DisplayName',[num2str(dias(k+18),'%.1f') ' \mum'])
	hold on;     
end
legend('NumColumns',2,'Location','West')
axis tight
ax=gca;
ax.XGrid='off';
ylabel({'IFCB';'[# L^-^1]'})
set(gca,'XTick',ticks,'XTickLabel',tickLabels,'FontSize',14);
xlim([datenum(2017,6,26,12,0,0) datenum(2017,7,14,12,0,0)]);
vline([t(cyc(1)) t(cyc(end)) t(acyc(1)) t(acyc(end))],'k:')
xlabel('Jun-Jul, 2017')
addlabel('d)')
set(gcf,'Position',[0.1120    0.0630    0.3578    0.8352])
cd(figsdir)
savefig(gcf,'Fig3_TimeSeriesAll.fig','compact')
print(gcf,'Fig3_TimeSeriesAll.png','-dpng','-r300')



[r, p]=corrcoef(ecobbp,n_ideal,'rows','complete'); r=r(1,2).^2 % r=0.57
[r, p]=corrcoef(cstarcp,n_ideal,'rows','complete') ; r=r(1,2).^2 %% r=0.83
[r, p]=corrcoef(psdsum_ifcb,n_ideal,'rows','complete'); r=r(1,2).^2 % % r=0.22
[r, p]=corrcoef(bbpratio,n_ideal,'rows','complete'); r=r(1,2).^2 % % r=-0.57


nanmean(n_ideal) % 1.0499
nanstd(n_ideal) % 0.0026

nanmean(n_ideal(cyc)) % 1.0480
nanstd(n_ideal(cyc)) % 0.0019

nanmean(n_ideal(acyc))% 1.0513
nanstd(n_ideal(acyc)) % 0.0021

% eddy to eddy mean difference in n
((0.051-0.048)/0.048 )*100 % 6%


%% Fig. 4: 
clearvars -except datadir mdir procdata figsdir cyc acyc
cd(mdir)

Get_Integrated_PSD_IFCB_perEddy %this is similar to what was done in Get_CruiseIntegrated_PSD_IFCB, but for the different eddies and for sunrise vs sunset.
ifcbPsdnFull_srise(ifcbPsdnFull_srise==0)=NaN;
ifcbPsdnFull_sset(ifcbPsdnFull_sset==0)=NaN;

load([procdata '\ifcb_IntegratedHourly.mat'],'ifcbPsdnBinHr','ifcbPsdBinHr','*Time*');


tzone = 10;
t=ifcbTimeHr-tzone/24;
tt=datevec(t);
y=ones(1,18).*2017; m = [6 6 6 6 7 7 7 7 7 7 7 7 7 7 7 7 7 7]; d=[27:30 1:14];
ticks = datenum(y,m,d); 
tickLabels = d;
linewidth = 1;



% Slope: 4.7394, SE: 0.1182, STD: 0.4094
% N0:10^4.9046, SE of exponent: 0.0552, STD of exponent: 0.1912
% 
% N0(+1sigma): 1.246765e+05
% N0: 8.028211e+04
% N0(-1sigma): 5.169553e+04
 
N0 = 8.028211e+04; % is it OK that this number is not the same as the average N0 at 4.3 um
Nminus = 5.169553e+04;
Nplus = 1.246765e+05;

N = N0.*((dias./dias(19)).^(-4.7));
N1 = N0.*((dias./dias(19)).^(-4.3));
N2 = N0.*((dias./dias(19)).^(-5.1));
N1b = Nminus*((dias./dias(19)).^(-4.7));
N2b = Nplus.*((dias./dias(19)).^(-4.7));


figure 
subplot 121
h1=loglog(dias(19:end),ifcbPsdnFull(19:end),'ko','MarkerFaceColor','w'); hold on
h2=loglog(dias(19:end),ifcbPsdnFull_cyc(19:end),'ko','MarkerFaceColor',rgb('blue')); hold on
h3=loglog(dias(19:end),ifcbPsdnFull_acyc(19:end),'ko','MarkerFaceColor',rgb('light red')); hold on
h4=loglog(dias(1:18),N(1:18),'k.');
h6=loglog(dias(1:18),N1(1:18),':k');
loglog(dias(1:18),N2(1:18),':k');
loglog(dias(1:18),N1b(1:18),':k');
loglog(dias(1:18),N2b(1:18),':k');
axis tight
h5=plot(0.6,(1.8*10^8)/0.5,'^k','MarkerFaceColor','g','MarkerSize',5) % pro+hetbackt 7*10^8 part/L to part/L/um divided by 0.6um if assuming that pro+hetbac average range 0.2-0.8um in size
ylabel('particles L^-^1 \mum^-^1')
axis tight
addlabel('a)');
grid on
xlim([0 100])
ylim([10^-1.6 10^12.5])
xlabel('Diameter [\mum]')
legend([h1(1) h2(1) h3(1) h4 h5],'cruise average','cyclonic eddy','anticyclonic eddy','extrapolation','{\it{Prochl.}} ref')
set(gca,'FontSize',14)
subplot 122
h1=loglog(dias(19:end),ifcbPsdnFull_srise(19:end),'ko','MarkerFaceColor','y'); hold on
h2=loglog(dias(19:end),ifcbPsdnFull_sset(19:end),'ko','MarkerFaceColor',rgb('purple')); hold on
axis tight
legend([h1(1) h2(1)],'7:00 hr','20:00 hr')
grid on
xlim([0 100])
ylim([10^-0.5 10^5.5])
xlabel('Diameter [\mum]')
ylabel('particles L^-^1 \mum^-^1')
addlabel('b)');

cd(figsdir)
savefig(gcf,'Fig4_PSD.fig','compact')
print(gcf,'Fig4_PSD.png','-dpng','-r300')


% Get log-log regression of full cruise PSD
indBin   = [19:30]; % 19 is the first good bin; bin 30 is the last bin with at least ~10 samples...bins>30 are likely still underestimating counts at those sizes
mdl      = fitlm(log10(dias(indBin)')-log10(4.3279),log10(ifcbPsdnFull_srise(indBin)'));
modelfun = @(b,x)b(1).*(x./4.3279).^b(2);
x1       = logspace(log10(0.2),log10(4.3279),50);
y1       = modelfun([10^mdl.Coefficients{1,1},mdl.Coefficients{2,1}],x1);
plot(x1,y1,'k-')
x1       = logspace(log10(3),log10(70),50);
y1       = modelfun([10^mdl.Coefficients{1,1},mdl.Coefficients{2,1}],x1);
figure;plot(log10(dias(indBin)'),log10(ifcbPsdnFull(indBin)'),'.')
hold on;plot(log10(x1),log10(y1),'.')

% Print out report
out = mdl.Coefficients{1:2,1:2};
fprintf('\nSlope: %0.4f, SE: %0.4f, STD: %0.4f\nN0:10^%0.4f, SE of exponent: %0.4f, STD of exponent: %0.4f\n',...
    -1*out(2),out(4),out(4)*sqrt(length(indBin)),out(1),out(3),out(3)*sqrt(length(indBin)))
fprintf('\nN0(+1sigma): %i\nN0: %i\nN0(-1sigma): %i\n\n',10^(out(1)+out(3)*sqrt(length(indBin))),...
    10^out(1),10^(out(1)-out(3)*sqrt(length(indBin))))


% Get log-log regression of full cruise PSD
indBin   = [19:30]; % 19 is the first good bin; bin 30 is the last bin with at least ~10 samples...bins>30 are likely still underestimating counts at those sizes
mdl      = fitlm(log10(dias(indBin)')-log10(4.3279),log10(ifcbPsdnFull_sset(indBin)'));
modelfun = @(b,x)b(1).*(x./4.3279).^b(2);
x1       = logspace(log10(0.2),log10(4.3279),50);
y1       = modelfun([10^mdl.Coefficients{1,1},mdl.Coefficients{2,1}],x1);
plot(x1,y1,'k-')
x1       = logspace(log10(3),log10(70),50);
y1       = modelfun([10^mdl.Coefficients{1,1},mdl.Coefficients{2,1}],x1);
figure;plot(log10(dias(indBin)'),log10(ifcbPsdnFull(indBin)'),'.')
hold on;plot(log10(x1),log10(y1),'.')

% Print out report
out = mdl.Coefficients{1:2,1:2};
fprintf('\nSlope: %0.4f, SE: %0.4f, STD: %0.4f\nN0:10^%0.4f, SE of exponent: %0.4f, STD of exponent: %0.4f\n',...
    -1*out(2),out(4),out(4)*sqrt(length(indBin)),out(1),out(3),out(3)*sqrt(length(indBin)))
fprintf('\nN0(+1sigma): %i\nN0: %i\nN0(-1sigma): %i\n\n',10^(out(1)+out(3)*sqrt(length(indBin))),...
    10^out(1),10^(out(1)-out(3)*sqrt(length(indBin))))


%% Fig. 4-5-6 (n plots and Sensitivity analysis plots)
plot_all_n

cd(procdata)



