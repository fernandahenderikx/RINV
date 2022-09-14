% plot all n_ideals 
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


cd(procdata)
filelist = dir('n_slope*'); % everytime you run additional analysis in MAIN.m the order of these may change. Be careful when calling the appropriate file below.
load('cyc_acyc_ind.mat'); % indices for start/end of sampling in the cyclone or anticyclone


% plot all at a glance, get the mean n and % diel change in n per filelist iteration
    figure_vertical
    ha = tight_subplot(4,3,[.05 .08],[.1 .2],[.2 .01])
for k = [1:length(filelist)]
    load(filelist(k).name)
axes(ha(k))
a = nan(18,24);
lt = gmt-10/24;
x = floor(lt);
[u uu]=unique(x);
clear a n_dailymean n_dailystd dayindex
for i = 1:length(u)
    ind = find(x==u(i)); 
    t = datevec(lt(ind));
    idx = t(:,4)+1; % to go from 1 to 24
    a(i,idx) = n_ideal(ind);
     
end
a(a==0)=NaN;
avgn{k,:} = char(sprintfc('%0.3f',round(nanmean(n_ideal),3))); 
stdevn{k,:} = char(sprintfc('%0.3f',round(nanstd(n_ideal),3))); 

deltaavgn{k,:} = char(sprintfc('%0.1f',round(nanmean((((a([8:10 13:15],20)-1)- (a([8:10 13:15],7)-1) )./ (a([8:10 13:15],7)-1)))*100,1))); % % 11.4, 9.7, 10.7, 9.5 
deltastdevn{k,:} =  char(sprintfc('%0.1f',round(nanstd((((a([8:10 13:15],20)-1)- (a([8:10 13:15],7)-1) )./ (a([8:10 13:15],7)-1)))*100,1))); % % 3.8, 2.4, 3.3, 2.6 

h1=plot([0.5:1:23.5],a([1:7 15:end],:)','ok','MarkerFaceColor',[.7 .7 .7],'DisplayName','outside eddies'); hold on
plot([0.5:1:23.5],a([11:12],:)','ok','MarkerFaceColor',[.7 .7 .7]); hold on
h2=plot([0.5:1:23.5],a([8:10],:)','ok','MarkerFaceColor','b','DisplayName','cyclone'); hold on
h3=plot([0.5:1:23.5],a([13:15],:)','ok','MarkerFaceColor','r','DisplayName','anticyclone'); hold on
title(filelist(k).name)
axis tight
xlim([0 24])
%ylim([1.01 1.085])
end



% load(filelist(3).name)
% clear M
% tzone = 10;
% M.t=gmt-tzone/24;
% tt=datevec(M.t);
% M.daytime = find(tt(:,4)>=7&tt(:,4)<=18);
% M.nighttime = find(tt(:,4)<=4|tt(:,4)>=19);
% y=ones(1,18).*2017; m = [6 6 6 6 7 7 7 7 7 7 7 7 7 7 7 7 7 7]; d=[27:30 1:14];
% M.ticks = datenum(y,m,d); 
% M.tickLabels = d;
% 
% lt = gmt-10/24;
% x = floor(lt);
% [u uu]=unique(x);
% clear a n_dailymean n_dailystd dayindex
% dayindex = nan(415,1);
% n_ideal_clean=n_ideal.*NaN;
% for i = 1:length(u)
%     ind = find(x==u(i)); if sum(~isnan(n_ideal(ind)))>16
%     n_dailymean(i,:) = nanmean(n_ideal(ind));
%     n_dailystd(i,:) = nanstd(n_ideal(ind));
%     n_ideal_clean(ind)=n_ideal(ind);
%     dayindex(ind) = i;
%     else
%         n_dailymean(i,[1:24])=NaN;
%         n_dailystd(i,[1:24])=NaN;
%         dayindex(ind) = i;
%         n_ideal_clean(ind)=NaN;          
%     end     
% end
% 
% 
% dm = nanmean(n_dailymean')-1
% ((nanmax(dm)-nanmin(dm))/nanmin(dm))*100 % 10%


%% Fig. 5 (results for slope = 4.57, No=fit)

k = [3]
    load(filelist(k).name)
    
clear a 

a = nan(18,24);
lt = gmt-10/24;
x = floor(lt);
[u uu]=unique(x);
for i = 1:length(u)
    ind = find(x==u(i)); 
    t = datevec(lt(ind));
    idx = t(:,4)+1; % to go from 1 to 24
    a(i,idx) = n_ideal(ind);
     
end

tt = lt; ttt = unique(floor(tt));
tstring =[datestr(ttt',6)];
% cyc = [346:537];
% acyc = [588:801]; % original
% transit = [1:345 538:587 801:926];
figure
subplot(2,2,[1:2])
plot(M.t, n_ideal,'ko-','MarkerSize',4); hold on
ylabel('Real Index of Refraction [{\it{n}}]')
ylim([1.022 1.032])
set(gca,'XTick',M.ticks,'XTickLabel',M.tickLabels,'FontSize',14);
xlim([datenum(2017,6,26,12,0,0) datenum(2017,7,14,12,0,0)]);
vline([M.t(cyc(1)) M.t(cyc(end)) M.t(acyc(1)) M.t(acyc(end))],'k--')
addlabel('a)')
text(M.t(cyc(1))+1,1.061,'cyclone')
text(M.t(acyc(1))+1,1.061,'anticyclone')
xlabel('Jun-Jul 2017')
subplot 224
% plot([1:24],a_percent([10 11],:)','ko:','MarkerFaceColor','b'); hold on
% plot([1:24],a_percent([15 16],:)','ko:','MarkerFaceColor','r'); hold on
h1=plot([0.5:1:23.5],a([1:7 15:end],:)','ok','MarkerFaceColor',[.7 .7 .7],'DisplayName','outside eddies'); hold on
plot([0.5:1:23.5],a([11:12],:)','ok','MarkerFaceColor',[.7 .7 .7]); hold on
h2=plot([0.5:1:23.5],a([8:10],:)','ok','MarkerFaceColor','b','DisplayName','cyclone'); hold on
h3=plot([0.5:1:23.5],a([13:15],:)','ok','MarkerFaceColor','r','DisplayName','anticyclone'); hold on
ylabel('Real Index of Refraction [{\it{n}}]')
xlim([0 24])
ylim([1.019 1.035])
text(2.8,1.034,['avg  = ' avgn{k,:} '{\pm}' stdevn{k,:} '%']);text(2.8,1.033,['{\Delta}  = ' deltaavgn{k,:} '{\pm}' deltastdevn{k,:} '%']);
legend([h1(1) h2(1) h3(1)],'Location','SouthEast')
addlabel('c)')
xlabel('Local Hour')
set(gca,'XTick',[3 6 9 12 15 18 21],'XTickLabel',{'3','6','9','12','15','18','21'})
subplot 223
boxplot(n_ideal_clean,dayindex,'plotstyle','compact','labels',tstring,'labelorientation','inline','colors',rgb('black')); hold on
ylim([1.022 1.032])
grid on
addlabel('b)')
ylabel('Real Index of Refraction [{\it{n}}]')
xlabel('Date [2017]')
set(gcf,'Position',[680   154   826   824])
cd(figsdir)
savefig(gcf,'Fig5_indexOfRefraction.fig','compact')
print(gcf,'Fig5_indexOfRefraction.png','-dpng','-r300')


% average daily dawn to dusk percent amplitude inside eddies
nanmean((((a([8:10 13:15],20)-1)- (a([8:10 13:15],7)-1) )./ (a([8:10 13:15],7)-1)))*100 % 11.3%
dm = n_dailymean-1;
((nanmax(dm)-nanmin(dm))/nanmin(dm))*100 % 10.1



%% Fig. 6 (results for slope = 3.91 & No=fit, slope = 5.23 & No=fit; slope = 4.57 & No fit - 1 stdev; slope = 4.57 & No fit + 1 stdev)

% fig 6 a, b, c, d
figure_vertical
    ha = tight_subplot(2,2,[.05 .08],[.1 .2],[.2 .01])
    f=1;
for k = [1 10 8 9] % make sure fig text(...) matches these numbers.
    load(filelist(k).name)
axes(ha(f))
a = nan(18,24);
lt = gmt-10/24;
x = floor(lt);
[u uu]=unique(x);
clear a n_dailymean n_dailystd dayindex
for i = 1:length(u)
    ind = find(x==u(i)); 
    t = datevec(lt(ind));
    idx = t(:,4)+1; % to go from 1 to 24
    a(i,idx) = n_ideal(ind);
     
end
a(a==0)=NaN;

h1=plot([0.5:1:23.5],a([1:7 15:end],:)','ok','MarkerFaceColor',[.7 .7 .7],'DisplayName','outside eddies'); hold on
plot([0.5:1:23.5],a([11:12],:)','ok','MarkerFaceColor',[.7 .7 .7]); hold on
h2=plot([0.5:1:23.5],a([8:10],:)','ok','MarkerFaceColor','b','DisplayName','cyclone'); hold on
h3=plot([0.5:1:23.5],a([13:15],:)','ok','MarkerFaceColor','r','DisplayName','anticyclone'); hold on
%title(filelist(k).name)
axis tight
xlim([0 24])
%ylim([1.01 1.085])
f=f+1;
end
axes(ha(1))
set(gca,'XTick',[3 6 9 12 15 18 21],'XTickLabel',''); 
legend([h1(1) h2(1) h3(1)])
ylim([1.028 1.060])
addlabel('a)')
text(3,1.058,'{\xi = 3.91}');
text(3,1.0555,'{N_o  = fit}');
text(3,1.0535,['{avg  = ' avgn{1,:} '{\pm}' stdevn{1,:} '%}']);
text(3,1.0515,['{\Delta  = ' deltaavgn{1,:} '{\pm}' deltastdevn{1,:} '%}']);
ylabel('Real Index of Refraction [{\it{n}}]')
axes(ha(2))
set(gca,'XTick',[3 6 9 12 15 18 21],'XTickLabel','')
%set(gca,'YTickLabel','')
ylim([1.01 1.022])
addlabel('b)')
text(3,1.0215,'{\xi = 5.23}'); text(3,1.0205,'{N_o  = fit}')
text(3,1.01955,['avg  = ' avgn{10,:} '{\pm}' stdevn{10,:} '%']);
text(3,1.01855,['{\Delta}  = ' deltaavgn{10,:} '{\pm}' deltastdevn{10,:} '%']);
axes(ha(3))
ylabel('Real Index of Refraction [{\it{n}}]')
set(gca,'XTick',[3 6 9 12 15 18 21],'XTickLabel',[3 6 9 12 15 18 21])
addlabel('c)')
xlabel('Local Hour');
ylim([1.02 1.045])
text(3,1.0435,'{\xi = 4.57}'); 
text(3,1.0415,'{N_o = fit - 1 \sigma}');
text(3,1.0395,['avg  = ' avgn{8,:} '{\pm}' stdevn{8,:} '%']);
text(3,1.0375,['{\Delta}  = ' deltaavgn{8,:} '{\pm}' deltastdevn{8,:} '%']);
axes(ha(4))
set(gca,'XTick',[3 6 9 12 15 18 21],'XTickLabel',[3 6 9 12 15 18 21])
%set(gca,'YTickLabel','')
addlabel('d)')
ylim([1.018 1.035])
text(3,1.034,'{\xi = 4.57}'); 
text(3,1.0325,'{N_o  = fit + 1 \sigma}')
text(3,1.031,['avg  = ' avgn{9,:} '{\pm}' stdevn{9,:} '%']);
text(3,1.0295,['{\Delta}  = ' deltaavgn{9,:} '{\pm}' deltastdevn{9,:} '%']);
xlabel('Local Hour');
cd(figsdir)
savefig(gcf,'Fig6_SensitivityAnalysis_1.fig','compact')
print(gcf,'Fig6_SensitivityAnalysis_1.png','-dpng','-r300')

%Fig. 6e; (results for slope = 4.57 but letting bins 4-6 um oscillate over the diel)
    figure
    k = [2]
    load(filelist(k).name)
a = nan(18,24);
lt = gmt-10/24;
x = floor(lt);
[u uu]=unique(x);
clear a_* c_* n_dailymean n_dailystd dayindex
for i = 1:length(u)
    ind = find(x==u(i)); 
    t = datevec(lt(ind));
    idx = t(:,4)+1; % to go from 1 to 24
    a(i,idx) = n_ideal(ind);
end
h1=plot([0.5:1:23.5],a([1:7 15:end],:)','ok','MarkerFaceColor',[.7 .7 .7],'DisplayName','outside eddies'); hold on
plot([0.5:1:23.5],a([11:12],:)','ok','MarkerFaceColor',[.7 .7 .7]); hold on
h2=plot([0.5:1:23.5],a([8:10],:)','ok','MarkerFaceColor','b','DisplayName','cyclone'); hold on
h3=plot([0.5:1:23.5],a([13:15],:)','ok','MarkerFaceColor','r','DisplayName','anticyclone'); hold on
%title(filelist(k).name)
axis tight
xlabel('Local Hour'); xlim([0 24])
ylim([1.018 1.04])
set(gca,'XTick',[3 6 9 12 15 18 21])
%legend([h1(1) h2(1) h3(1)])
text(3,1.039,'{\xi = 4.57}'); text(3,1.0375,'{N_o = variable}')
text(3,1.036,['avg  = ' avgn{k,:} '{\pm}' stdevn{k,:} '%']);text(3,1.0345,['{\Delta}  = ' deltaavgn{k,:} '{\pm}' deltastdevn{k,:} '%']);
addlabel('e)')
ylabel('Real Index of Refraction [{\it{n}}]')
cd(figsdir)
set(gcf,'Position',[2815         323         280         366])
savefig(gcf,'Fig6_SensitivityAnalysis_2.fig','compact')
print(gcf,'Fig6_SensitivityAnalysis_2.png','-dpng','-r300')

% average daily dawn to dusk percent amplitude inside eddies
nanmean((((a([8:10 13:15],20)-1)- (a([8:10 13:15],7)-1) )./ (a([8:10 13:15],7)-1)))*100 % 10.1%


%% sensitivity analysis - truncate data at 10, 20, 50 um
  figure_vertical
    ha = tight_subplot(1,3,[.05 .08],[.1 .2],[.2 .01])
    f=1;
for k = [4 5 6] % make sure these numbers match the numbers called in text(...) below
    load(filelist(k).name)
axes(ha(f))
clear a2 n_dailymean n_dailystd dayindex
a2 = nan(18,24);
lt = gmt-10/24;
x = floor(lt);
[u uu]=unique(x);

for i = 1:length(u)
    ind = find(x==u(i)); 
    t = datevec(lt(ind));
    idx = t(:,4)+1; % to go from 1 to 24
    a2(i,idx) = n_ideal(ind);
     
end

h1=plot([0.5:1:23.5],a([1:7 15:end],:)','ok','MarkerFaceColor',[.7 .7 .7],'DisplayName','outside eddies'); hold on
plot([0.5:1:23.5],a([11:12],:)','ok','MarkerFaceColor',[.7 .7 .7]); hold on
h2=plot([0.5:1:23.5],a([8:10],:)','ok','MarkerFaceColor','b','DisplayName','cyclone'); hold on
h3=plot([0.5:1:23.5],a([13:15],:)','ok','MarkerFaceColor','r','DisplayName','anticyclone'); hold on
%title(filelist(k).name)
 xlim([0 24])
ylim([1.02 1.035])
f=f+1;
end
axes(ha(1))
set(gca,'XTick',[3 6 9 12 15 18 21],'XTickLabel',''); 
legend([h1(1) h2(1) h3(1)])
addlabel('a)')
text(3,1.034,'{\xi = 4.57, max 10 \mum}');
text(3,1.033,['{avg  = ' avgn{4,:} '{\pm}' stdevn{4,:} '%}']);
text(3,1.032,['{\Delta  = ' deltaavgn{4,:} '{\pm}' deltastdevn{4,:} '%}']);
ylabel('Real Index of Refraction [{\it{n}}]')
axes(ha(2))
set(gca,'XTick',[3 6 9 12 15 18 21],'XTickLabel','')
set(gca,'YTickLabel','')
addlabel('b)')
text(3,1.034,'{\xi = 4.57, max = 20 \mum}'); text(3,1.033,'{N_o  = fit}')
text(3,1.032,['avg  = ' avgn{5,:} '{\pm}' stdevn{5,:} '%']);
text(3,1.031,['{\Delta}  = ' deltaavgn{5,:} '{\pm}' deltastdevn{5,:} '%']);
axes(ha(3))
ylabel('Real Index of Refraction [{\it{n}}]')
set(gca,'XTick',[3 6 9 12 15 18 21],'XTickLabel',[3 6 9 12 15 18 21])
addlabel('c)')
xlabel('Local Hour');
text(3,1.034,'{\xi = 4.57, max = 50 \mum}'); 
text(3,1.033,['avg  = ' avgn{6,:} '{\pm}' stdevn{6,:} '%']);
text(3,1.032,['{\Delta}  = ' deltaavgn{6,:} '{\pm}' deltastdevn{6,:} '%']);
cd(figsdir)
savefig(gcf,'Fig6_SensitivityAnalysis_3.fig','compact')
print(gcf,'Fig6_SensitivityAnalysis_3.png','-dpng','-r300')





