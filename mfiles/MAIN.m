%% MAIN script for Index of Refraction paper. Certain analysis that take a while are commented out below, but this mfile reproduces all analysis and all figures from the manuscript.

% main directory RINV_v3 should contain the following folders: mfiles, data, proc, figures.

clear all; close all;
mdir=cd;
cd ..\..
%addpath(genpath('RINV_v3'))
cd(mdir)
cd ..
cd data
datadir = cd;
cd ..\proc
procdata = cd;
cd ..\figures
figsdir = cd;
cd(mdir)


%% process IFCB data (aggregated cruise PSD). Get parameters of power law fit for extrapolation to 0.2 um.
cd(mdir)
Get_CruiseIntegrated_PSD_IFCB % generates ifcbPSDfitReport.txt report; also report generated at the end of script which shows slope and No for uncertainty analysis of the extrapolated portion of the PSD

% result pasted here for ease:
% assuming no weighting function
% Slope: 4.7907, SE: 0.0749, STD: 0.2370
% N0:10^4.7481, SE of exponent: 0.0286, STD of exponent: 0.0903
% 
% N0(+1sigma): 6.892796e+04
% N0: 5.598970e+04
% N0(-1sigma): 4.548004e+04

% assuming weights = number of particles
% Slope: 4.5694, SE: 0.2098, STD: 0.6635
% N0:10^4.7252, SE of exponent: 0.0288, STD of exponent: 0.0909
% 
% N0(+1sigma): 6.547742e+04
% N0: 5.310886e+04
% N0(-1sigma): 4.307670e+04
% 
% 
cd(mdir)
Get_HourlyIntegrated_PSD_IFCB % get 415x38 PSDs. This generates the main IFCB-based dataset used in the paper (no extrapolation to 0.2 um yet).

%% load and process min avgd beamc, cleanup outliers if any, bin hourly, make size match IFCB hourly.
cd(mdir)
CleanupIOPs % get binData with hourly time, cp, bbp

%% Get reference diameter per bin for the average slope and No. refdiam is recalculated within each of the runs below when a different slope or No is used.
get_refdiam

%%  Estimate n based on different assumptions for the extrapolated portion of the IFCB
cd(mdir)

% fixed parameters
ifcbdata = [procdata '\ifcb_IntegratedHourly.mat'];
load([procdata '\cp_bbp_hourly.mat'],'cstarcp');
load([procdata '\refdiam.mat'],'refdiam');  % here only really care for the bin number of the ref, not the mean values per bin. Those will be recalculated within each run to account for the different slopes and No.

%refwater = get_water_n(660,26,35) % == 1.337
refwater=1.337; % Water refractive index
lambda = 660; % nm; will convert to meter inside code below
% % choice of n':
% lambdaair=lambda/1.337; % in water
% nimag= 0.010658*exp(-0.007186*lambda) %  == 0.0003; Stramski et al 2004
nimag = 0.0003;

tic % this takes a long time. If need to speed it up, especially for testing purposes, change increments of nall and nint in calc_cp.m
% % assume everything below refdiam is static in time.
slope = 4.57;
No = 10^4.7252;
lutfile = [procdata '\LUTcp_slope457_NoAvg.mat']; % name to save as
nfile = [procdata '\n_slope457_NoAvg.mat']; % name to save as
run2a(ifcbdata,slope,No,refwater,nimag,lambda,refdiam,cstarcp,lutfile,nfile)

% change slopes
slope = 3.91;
No = 10^4.7252;
lutfile = [procdata '\LUTcp_slope391_NoAvg.mat']; % name to save as
nfile = [procdata '\n_slope391_NoAvg.mat']; % name to save as
run2a(ifcbdata,slope,No,refwater,nimag,lambda,refdiam,cstarcp,lutfile,nfile)

slope = 5.23;
No = 10^4.7252;
lutfile = [procdata '\LUTcp_slope523_NoAvg.mat']; % name to save as
nfile = [procdata '\n_slope523_NoAvg.mat']; % name to save as
run2a(ifcbdata,slope,No,refwater,nimag,lambda,refdiam,cstarcp,lutfile,nfile)

% change No
slope = 4.57;
No = 4.307670e+04;
lutfile = [procdata '\LUTcp_slope457_NoMinusOneSigma.mat']; % name to save as
nfile = [procdata '\n_slope457_NoMinusOneSigma.mat']; % name to save as
run2a(ifcbdata,slope,No,refwater,nimag,lambda,refdiam,cstarcp,lutfile,nfile)

slope = 4.57;
No = 6.547742e+04;
lutfile = [procdata '\LUTcp_slope457_NoPlusOneSigma.mat']; % name to save as
nfile = [procdata '\n_slope457_NoPlusOneSigma.mat']; % name to save as
run2a(ifcbdata,slope,No,refwater,nimag,lambda,refdiam,cstarcp,lutfile,nfile)

% % assume 4:6 oscillate WITH refdiam. everything else static.
slope = 4.57;
No = 10^4.7252;
lutfile = [procdata '\LUTcp_slope457_4to6diel_NoAvgmat']; % name to save as
nfile = [procdata '\n_slope457_4to6diel_NoAvg.mat']; % name to save as
run2b(ifcbdata,slope,No,refwater,nimag,lambda,refdiam,cstarcp,lutfile,nfile)
toc

% % additional sensitivity analysis 
% truncate main analysis to 10um; 20 um; 50 um; 6 um to confirm that large
% particles do not contribute with n estimates (both because of how Qc
% varies with diameter for differnet n in the > 10 um size region, and also
% because particles > 10 um make up a negligible portion of the cp signal
% to begin with.
tic
slope = 4.57;
No = 10^4.7252;
lutfile = [procdata '\LUTcp_slope457_NoAvg.mat']; % nfile name will change within run3 
nfile = [procdata '\n_slope457_NoAvg.mat']; % nfile name will change within run3 
run3(ifcbdata,slope,No,refwater,nimag,lambda,refdiam,cstarcp,lutfile,nfile)
toc

% % assume 2:6 oscillate WITH 6um. everything else static.
% slope = 4.57;
% No = 10^4.7252;
% lutfile = [procdata '\LUTcp_slope457_2to6diel_NoAvg.mat']; % name to save as
% nfile = [procdata '\n_slope457_2to6diel_NoAvg.mat']; % name to save as
% run2c(ifcbdata,slope,No,refwater,nimag,lambda,refdiam,cstarcp,lutfile,nfile)




%% APPENDIX analysis

% %1) get LISST pscat
% %create LISST_pscat.mat. 

cd(datadir)
cd LISST\mat_files
load('MINAVG_LISST_CORR.mat') % min averaged LISST data corrected for filtered sections (see Henderikx-Freitas et al 2020, Applied Optics)
cd ..\mfiles
get_pscat(all_min_time,all_min_data) % make sure this uses James's modified pscat code (to output raw and tau)

%2) calculate Mie VSF for the PSD used in the paper 
%This creates vsf_ifc_Full_ALLn.mat. 
cd(mdir)
Get_MieVSF_forPaperPSD % takes ~80 seconds.
cd(procdata)
load vsf_ifcb_Full_ALLn.mat
load('LISST_pscat.mat','sdate','pscat_all')
pscat_all(pscat_all<=0)=NaN;

% determine a fudge factor to convert LISST pscat to VSF. This is
% arbitrary, and was chosen to allow retrieved LISST VSFs to generally match
% magnitude of Mie VSFs.
fudge_factor = 9.*10^-6; %coefficient [m-1 sr-1 count-1] relating the LISST measured pscat to mie vsf

% LISST pscat is in 1 min intervals. Make hourly averages here.
xi = gmt;
clear pscat_hr
for i = 1:length(xi)
    ind = find(sdate>=xi(i)-0.5/24&sdate<xi(i)+0.5/24);
    pscat_hr(i,:) = nanmean(pscat_all(ind,:));
end


% find minimum square error between LISST pscat*fudge factor (aka "VSF")
% and Mie VSF for the selected angles we have generated the data for
ang_target_idx = [14 18 20 21 22 23 24 25 26]; % these are the angle bins that we have generated the Mie VSF for
k=1;
figure
for i = 1:length(ang_target_idx)
    pscat_hr_corr = pscat_hr(:,ang_target_idx(i)).*fudge_factor;
    pscat_hr_corr=pscat_hr_corr';

    vsf_ifcb = eval(['vsf_ifcb_ang' num2str(ang_target_idx(i))]);

    mse = (repmat(pscat_hr_corr,length(nint),1) - vsf_ifcb).^2;
    
    subplot(3,3,i)
    plot(mse)
    title(['MSE for angle = ' num2str(ang_target(ang_target_idx(i)))]);
    
    [ind indd] = min(mse);
    nint(indd) % 1.0310
    
    lisst_nideal(k,:) = nint(indd);

    k=k+1;
end


% interp lisst n result to match the 415 from beamc etc.
cd(procdata)
D=load('n_slope457_NoAvg.mat')
for i = 1:length(ang_target_idx)    
    lisst_nideal_i(i,:) = interp1(gmt,lisst_nideal(i,:),D.gmt);
end

% put back nans. there's a better way to do this.
ind = isnan(D.n_ideal);
lisst_nideal_i(:,ind)=NaN;
% cleanup last lisst n just to make plotting easier - looks like an outlier
lisst_nideal_i(:,end)=NaN;

% choose an angle for plots/analysis
ang_target_single = 9; % 9 = angle 26.
vsf_ifcb = vsf_ifcb_ang26;

gmt_i = D.gmt;
save n_slope457_NoAvg_LISST.mat lisst_nideal_i ang_target_idx ang_target gmt_i


% Compare Mie VSF at selected sunrise/sunset example points with LISST
% pscat. pick a couple of n for the Mie VSF.
cd(mdir)
get_theoreticalVSF_sunrise_sunset
load vsf_ifcb_sunrise_sunset.mat % created in get_theoreticalVSF_sunrise_sunset

% make sure gmt below is the trimmed version (size 349 and not 415).
ind173 = find(sdate>=gmt(173)-0.5/24&sdate<=gmt(173)+0.5/24);
ind185 = find(sdate>=gmt(185)-0.5/24&sdate<=gmt(185)+0.5/24);

% Appendix Fig A1
figure
subplot 121
plot(ang_target,vsf_ifcb_sunrise_sunset_n102(:,1),'ob','MarkerSize',3,'MarkerFaceColor','b'); hold on
plot(ang_target,vsf_ifcb_sunrise_sunset_n102(:,2),'or','MarkerSize',3,'MarkerFaceColor','r')
plot(ang_target,vsf_ifcb_sunrise_sunset_n104(:,1),'ob','MarkerSize',3,'MarkerFaceColor','w'); 
plot(ang_target,vsf_ifcb_sunrise_sunset_n104(:,2),'or','MarkerSize',3,'MarkerFaceColor','w')
set(gca,'YScale','Log')
set(gca,'XScale','Log')
ylim([5*10^-3 10^2])
xlim([0 16])
legend('example sunrise PSD; {\it{n}} = 1.010', 'example sunset PSD; {\it{n}} = 1.010', 'example sunrise PSD; {\it{n}} = 1.025', 'example sunset PSD; {\it{n}} = 1.025','Location','SouthWest')
xlabel('Angle [degrees]')
ylabel('Theoretical VSF [m^-^1 sr^-^1]')
%vline(ang_target(26),'k--')
addlabel('a)')
set(gcf,'Position', [194         325        1138         580])
subplot 122
errorbar(ang_target,nanmean(pscat_all(ind173,:)),nanstd(pscat_all(ind173,:)),'ok','MarkerSize',3,'MarkerFaceColor',rgb('light blue')); hold on
errorbar(ang_target,nanmean(pscat_all(ind185,:)),nanstd(pscat_all(ind185,:)),'ok','MarkerSize',3,'MarkerFaceColor',rgb('light red')); hold on
set(gca,'YScale','Log')
set(gca,'XScale','Log')
ylim([5*10^3 10^8])
xlim([0 16])
legend('example sunrise PSD {\pm} 30 min','example sunset PSD {\pm} 30 min', 'Location','SouthWest')
xlabel('Angle [degrees]')
ylabel('LISST pscat [counts]')
vline(ang_target(26),'k--')
addlabel('b)')
cd(figsdir)
savefig(gcf,'FigSuppInfo1.fig','compact')
print(gcf,'FigSuppInfo1.png','-dpng','-r300')

% Appendix Fig A2
load('cyc_acyc_ind.mat')
clear M
tzone = 10;
M.t=D.gmt-tzone/24;
tt=datevec(M.t);
M.daytime = find(tt(:,4)>=7&tt(:,4)<=18);
M.nighttime = find(tt(:,4)<=4|tt(:,4)>=19);
y=ones(1,18).*2017; m = [6 6 6 6 7 7 7 7 7 7 7 7 7 7 7 7 7 7]; d=[27:30 1:14];
M.ticks = datenum(y,m,d); 
M.tickLabels = d;

figure_long
%subplot(1,4,[1:3])
subplot(1,3,[1:3])
plot(D.gmt-10/24,D.n_ideal,'ko-','MarkerSize',4,'MarkerFaceColor','k'); hold on
plot(D.gmt-10/24,lisst_nideal_i(ang_target_single,:),'ko-','MarkerSize',4,'MarkerFaceColor',rgb('cyan'))
ylim([1.018 1.033])
set(gca,'XTick',M.ticks,'XTickLabel',M.tickLabels,'FontSize',14);
xlim([datenum(2017,6,26,12,0,0) datenum(2017,7,14,12,0,0)]);
vline([M.t(cyc(1)) M.t(cyc(end)) M.t(acyc(1)) M.t(acyc(end))],'k--')
%addlabel('a)')
ylabel('Estimated real index of refraction [{\it{n}}]')
text(M.t(cyc(1))+1,1.0335,'cyclone')
text(M.t(acyc(1))+1,1.0335,'anticyclone')
xlabel('Jun-Jul 2017')
legend('c_p + IFCB estimated \it{n}','LISST + IFCB estimated \it{n}', 'Location','SouthEast')
cd(figsdir)
savefig(gcf,'FigSuppInfo2.fig','compact')
print(gcf,'FigSuppInfo2.png','-dpng','-r300')


% check sizes affecting angle 26.
cd(mdir)
check_sizes_affecting_angle26


%% Other analysis

% check possibility of multiple solutions when solving for n_ideal
 cd(procdata)
 load([procdata '\LUTcp_slope457_NoAvg.mat'],'cp','nrange','dias');
 load([procdata '\cp_bbp_hourly.mat'],'cstarcp');

[n_ideal sq] = Get_Estimated_n(cstarcp,cp,nrange);

figure_vertical
plot(nrange,sq,'.'); xlabel('index of refraction'); ylabel('least square difference')
ylim([-0.0005 0.008]); xlim([1 1.05])
title('slope = 4.57, No = avg; only one global minimum')

% PSD to decide minimum good bin
load([procdata '\ifcb_cruiseIntegrated.mat'],'ifcbPsdnFull');
figure 
loglog(dias(1:end),ifcbPsdnFull(1:end),'ko','MarkerFaceColor','w'); hold on
grid on
xlim([0 100])
ylim([10^-1.6 10^8])




%% FIGURES

close all; 

%% Fig. 1
cd(procdata)
clearvars -except datadir mdir procdata figsdir cyc acyc

load('LUTcp_slope457_NoAvg.mat')
load('refdiam.mat')

ind0 = find(round(nrange,4)==1.0200);
ind1 = find(round(nrange,4)==1.0500);
ind2 = find(round(nrange,4)==1.0800);
% 
% make up a PSD of slope 4
slope=4;
refdias = 21;
No = 46526.6616249321; 
N = No.*((dias(1:38)./dias(refdias)).^-slope);
N_m3 = N.*1000;
dias_m = dias.*10^-6; % meters
X0 = (pi/4.*(dias_m.^2)).*Qc(ind0,:).*N_m3;
X1 = (pi/4.*(dias_m.^2)).*Qc(ind1,:).*N_m3;
X2 = (pi/4.*(dias_m.^2)).*Qc(ind2,:).*N_m3;
X0_norm = X0./nansum(X0);
X1_norm = X1./nansum(X1);
X2_norm = X2./nansum(X2);
X0star = (pi/4.*(dias_m.^2)).*Qcstar(ind0,:).*N_m3;
X1star = (pi/4.*(dias_m.^2)).*Qcstar(ind1,:).*N_m3;
X2star = (pi/4.*(dias_m.^2)).*Qcstar(ind2,:).*N_m3;
X0star_norm = X0star./nansum(X0star);
X1star_norm = X1star./nansum(X1star);
X2star_norm = X2star./nansum(X2star);

figure
subplot 131
loglog(dias,Qc(ind0,:),'ko-','MarkerFaceColor','w','MarkerSize',5); hold on; 
loglog(dias,Qcstar(ind0,:),'k^-','MarkerFaceColor','w','MarkerSize',5); hold on; 
loglog(dias,Qc(ind1,:),'ko-','MarkerFaceColor','y','MarkerSize',5); hold on; 
loglog(dias,Qcstar(ind1,:),'k^-','MarkerFaceColor','y','MarkerSize',5); hold on; 
loglog(dias,Qc(ind2,:),'ko-','MarkerFaceColor','c','MarkerSize',5); hold on; 
loglog(dias,Qcstar(ind2,:),'k^-','MarkerFaceColor','c','MarkerSize',5); hold on; 
xlim([10^-1 10^2])
%vline(5.54,':k')
addlabel('a)')
ylabel('{\itQ_c or Q_c^*}   [unitless]')
legend('{\itQ_c; n} = 1.02','{\itQ_c^*; n} = 1.02','{\itQ_c; n} = 1.05','{\itQ_c^*; n} = 1.05','{\itQ_c; n} = 1.08','{\itQ_c^*; n} = 1.08','Location','SouthEast');
xlabel('Diameter [\mum]'); 
subplot 132
semilogx(dias,X0,'ko-','MarkerFaceColor','w','MarkerSize',5); hold on; 
semilogx(dias,X1,'ko-','MarkerFaceColor','y','MarkerSize',5); hold on; 
semilogx(dias,X2,'ko-','MarkerFaceColor','c','MarkerSize',5); hold on; 
semilogx(dias,X0star,'k^-','MarkerFaceColor','w','MarkerSize',5); hold on; 
semilogx(dias,X1star,'k^-','MarkerFaceColor','y','MarkerSize',5); hold on; 
semilogx(dias,X2star,'k^-','MarkerFaceColor','c','MarkerSize',5); hold on; 
axis tight
xlim([10^-1 10^2]); 
%vline(5.54,':k')
ylabel('c_p per size class [m^-^1]')
xlabel('Diameter [\mum]'); 
%hline(0,'k')
addlabel('b)')
subplot 133
semilogx(dias,cumsum(X0_norm),'ko-','MarkerFaceColor','w','MarkerSize',5); hold on; 
semilogx(dias,cumsum(X1_norm),'ko-','MarkerFaceColor','y','MarkerSize',5); hold on; 
semilogx(dias,cumsum(X2_norm),'ko-','MarkerFaceColor','c','MarkerSize',5); hold on; 
semilogx(dias,cumsum(X0star_norm),'k^-','MarkerFaceColor','w','MarkerSize',5); hold on; 
semilogx(dias,cumsum(X1star_norm),'k^-','MarkerFaceColor','y','MarkerSize',5); hold on; 
semilogx(dias,cumsum(X2star_norm),'k^-','MarkerFaceColor','c','MarkerSize',5); hold on; 
axis tight
xlim([10^-1 10^2]); 
ylim([0 1])
ylabel('norm. cum. dist. of c_p per size class')
xlabel('Diameter [\mum]'); 
addlabel('c)')
cd(figsdir)
set(gcf,'Position',[247         491        1317         411])
print(gcf,'Fig1_Qc.png','-dpng','-r300');
savefig(gcf,'Fig1_Qc.fig','compact')

% 
% figure; 
% semilogx(dias,cumsum(X0star_norm)); hold on % particles > 10 make up <5% of total cp for n=1.02
% semilogx(dias,cumsum(X1star_norm)); hold on % particles > 10 make up <1% of total cp for n=1.05
% semilogx(dias,cumsum(X2star_norm)); hold on % particles > 10 make up <1% of total cp for n=1.08
% 

% figure
% imagesc(log10(dias),nrange',Qc)
% colormap(cmocean('deep'))
% h=colorbar; title(h,'Q_c [unitless]')
% xlabel('Diameter [\mum]'); ylabel('real index of refraction [n]')
% set(gca,'XTick',[-0.5 0 0.5 1 1.5 2],'XTickLabel',{'0.3','1','3','10','30','100'},'FontSize',14);

%% Fig. 2: flowchart

%% Fig. 3
clearvars -except datadir mdir procdata figsdir cyc acyc
cd(procdata)
Get_Integrated_PSD_IFCB_perEddy %this is similar to what was done in Get_CruiseIntegrated_PSD_IFCB, but for the different eddies.
load([procdata '\ifcb_IntegratedHourly.mat'],'ifcbPsd*','*Time*','*Agg*');
load([procdata '\cp_bbp_hourly.mat']);
load([procdata '\LUTcp_slope457_NoAvg.mat'],'dias');
load([procdata '\n_slope457_NoAvg.mat']);

% indices for cyclone and anticyclone eddies:
cyc = [155:245]; 
acyc = [269:361];

psdsum_ifcb = nansum(ifcbPsdBinHr(:,21:38)');psdsum_ifcb(psdsum_ifcb==0)=NaN;

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
ha = tight_subplot(3,1,[.045 .03],[.08 .08],[.2 .2])
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
for k = 1 : 6
    plot(t,ifcbPsdnBinHr(:,k+20),'-', 'Color', listOfThermalColors(k, :), 'LineWidth', 1.5,'DisplayName',[num2str(dias(k+20),'%.2f') ' \mum'])
	hold on;     
end
legend('NumColumns',2,'Location','West')
axis tight
ax=gca;
ax.XGrid='off';
ylabel({'IFCB';'[# L^-^1 \mum^-^1]'})
set(gca,'XTick',ticks,'XTickLabel',tickLabels,'FontSize',14);
xlim([datenum(2017,6,26,12,0,0) datenum(2017,7,14,12,0,0)]);
vline([t(cyc(1)) t(cyc(end)) t(acyc(1)) t(acyc(end))],'k:')
xlabel('Jun-Jul, 2017')
xtickangle(0)
addlabel('c)')
set(gcf,'Position',[0.3312    0.1185    0.3906    0.7861])
cd(figsdir)
savefig(gcf,'Fig3_TimeSeriesAll.fig','compact')
print(gcf,'Fig3_TimeSeriesAll.png','-dpng','-r300')



[r, p]=corrcoef(cstarcp,n_ideal,'rows','complete') ; r=r(1,2).^2 %% r=0.75
[r, p]=corrcoef(psdsum_ifcb,n_ideal,'rows','complete'); r=r(1,2).^2 % % r=0.0


nanmean(n_ideal) % 1.0264
nanstd(n_ideal) % 0.0014

nanmean(n_ideal(cyc)) % 1.0254
nanstd(n_ideal(cyc)) % 0.0009

nanmean(n_ideal(acyc))% 1.0271
nanstd(n_ideal(acyc)) % 0.0010

% eddy to eddy mean difference in n
(((nanmean(n_ideal(acyc))-1) - (nanmean(n_ideal(cyc))-1)))./ (nanmean(n_ideal(cyc))-1)*100 % 6.8%


%% Fig. 4: 
clearvars -except datadir mdir procdata figsdir cyc acyc
cd(mdir)

Get_Integrated_PSD_IFCB_perEddy %this is similar to what was done in Get_CruiseIntegrated_PSD_IFCB, but for the different eddies and for sunrise vs sunset.
ifcbPsdnFull_srise(ifcbPsdnFull_srise==0)=NaN;
ifcbPsdnFull_sset(ifcbPsdnFull_sset==0)=NaN;
load([procdata '\ifcb_IntegratedHourly.mat'],'ifcbPsdnBinHr','ifcbPsdBinHr','*Time*');
load([procdata '\LUTcp_slope457_NoAvg.mat'],'dias');


tzone = 10;
t=ifcbTimeHr-tzone/24;
tt=datevec(t);
y=ones(1,18).*2017; m = [6 6 6 6 7 7 7 7 7 7 7 7 7 7 7 7 7 7]; d=[27:30 1:14];
ticks = datenum(y,m,d); 
tickLabels = d;
linewidth = 1;



% Slope: 4.5694, SE: 0.2098, STD: 0.6635
% N0:10^4.7252, SE of exponent: 0.0288, STD of exponent: 0.0909
% 
% N0(+1sigma): 6.547742e+04
% N0: 5.310886e+04
% N0(-1sigma): 4.307670e+04
 
N0 = 5.310886e+04; % is it OK that this number is not the same as the average N0 at 4.3 um
Nminus = 4.307670e+04;
Nplus = 6.547742e+04;

N = N0.*((dias./dias(21)).^(-4.57));
N1 = N0.*((dias./dias(21)).^(-3.91));
N2 = N0.*((dias./dias(21)).^(-5.23));
N1b = Nminus*((dias./dias(21)).^(-4.57));
N2b = Nplus.*((dias./dias(21)).^(-4.57));

refdias = 21;

figure 
subplot 121
h1=loglog(dias(refdias:end),ifcbPsdnFull(refdias:end),'ko','MarkerFaceColor','w','MarkerSize',5); hold on
h2=loglog(dias(refdias:end),ifcbPsdnFull_cyc(refdias:end),'ko','MarkerFaceColor',rgb('blue'),'MarkerSize',5); hold on
h3=loglog(dias(refdias:end),ifcbPsdnFull_acyc(refdias:end),'ko','MarkerFaceColor',rgb('light red'),'MarkerSize',5); hold on
h4=loglog(dias(1:refdias-1),N(1:refdias-1),'k.');
h6=loglog(dias(1:refdias-1),N1(1:refdias-1),':k');
loglog(dias(1:refdias-1),N2(1:refdias-1),':k');
loglog(dias(1:refdias-1),N1b(1:refdias-1),':k');
loglog(dias(1:refdias-1),N2b(1:refdias-1),':k');
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
h1=loglog(dias(refdias:end),ifcbPsdnFull_srise(refdias:end),'ko','MarkerFaceColor','y','MarkerSize',7); hold on
h2=loglog(dias(refdias:end),ifcbPsdnFull_sset(refdias:end),'ks','MarkerFaceColor',rgb('purple'),'MarkerSize',5); hold on
axis tight
legend([h1(1) h2(1)],'7:00 hr','20:00 hr')
grid on
xlim([0 100])
ylim([10^-0.5 10^5.5])
xlabel('Diameter [\mum]')
ylabel('particles L^-^1 \mum^-^1')
addlabel('b)');
set(gca,'FontSize',14)
set(gcf,'Position',[680   465   831   513]);
cd(figsdir)
savefig(gcf,'Fig4_PSD.fig','compact')
print(gcf,'Fig4_PSD.png','-dpng','-r300')

% srise slope: 4.58+-0.79
% sset slope: 4.56+-0.97


% Get log-log regression of srise psd
indBin   = [refdias:30]; % 19 is the first good bin; bin 30 is the last bin with at least ~10 samples...bins>30 are likely still underestimating counts at those sizes
mdl      = fitlm(log10(dias(indBin)')-log10(dias(refdias)),log10(ifcbPsdnFull_srise(indBin)'),'Weights',ifcbCountFull_srise(indBin));
modelfun = @(b,x)b(1).*(x./dias(refdias)).^b(2);
x1       = logspace(log10(0.2),log10(dias(refdias)),50);
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


% Get log-log regression of sset psd
indBin   = [refdias:30]; % 19 is the first good bin; bin 30 is the last bin with at least ~10 samples...bins>30 are likely still underestimating counts at those sizes
mdl      = fitlm(log10(dias(indBin)')-log10(dias(refdias)),log10(ifcbPsdnFull_sset(indBin)'),'Weights',ifcbCountFull_sset(indBin));
modelfun = @(b,x)b(1).*(x./dias(refdias)).^b(2);
x1       = logspace(log10(0.2),log10(dias(refdias)),50);
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




%% Fig. 5-6 (n plots and Sensitivity analysis plots)
plot_all_n

cd(mdir)



