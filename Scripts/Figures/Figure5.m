
addpath(genpath(pwd));
%addpath(genpath('/Volumes/MELANOMAII'));

saveFigures = true;

%% Figure5D

clearvars -except saveFigures
clc
% 
% load('/Volumes/MELANOMA/Data/Simulations/5nodes/S_outpar1000_5_1_92')
% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Simulations/5nodes/S_outpar1000_5_1_92')
load('./Data/Data1000')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;
param = 5;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
% rectangle('Position',[606150 0 170 300],'FaceColor',[0.8 .8 .8],'EdgeColor','none')
% hold on
h = plot(969400:970400,S_outpar{param}(1:5,969400:970400),'LineWidth',2);
hold on 
plot(969400:970400,repmat(thres(915),1,length(969400:970400)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'}, {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255;...
    [44 90 160]./255; [33 68 120]./255});
xlim([969400, 970400])
xticks([969400, 970400])
xticklabels({'0', '1000'})
yticks([0, 250])
yticklabels({'0', '2.5'})
ylim([0,250])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
ylabel('gene product (10^2)')
xlabel('time')
box off

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure5/Figure5D'])
    print('-dpdf',['./Figures/Figure5/Figure5D'])
end

%% Figure5E

clearvars -except saveFigures
clc

% load('/Volumes/MELANOMA/Data/Simulations/5nodes/S_outpar1000_5_10_92')
% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Simulations/5nodes/S_outpar1000_5_10_92')
load('./Data/Data1000')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;
param = 5;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(41000:42000,S_outpar{param}(1:5,41000:42000),'LineWidth',2);
hold on 
plot(41000:42000,repmat(thres(915),1,length(41000:42000)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'}, {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255;...
    [44 90 160]./255; [33 68 120]./255});
xlim([41000, 42000])
xticks([41000, 42000])
xticklabels({'0', '1000'})
yticks([0, 250])
yticklabels({'0', '2.5'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
ylabel('gene product (10^2)')
xlabel('time')
ylim([0,250])
box off

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure5/Figure5E'])
    print('-dpdf',['./Figures/Figure5/Figure5E'])
end

%% Figure5F

clearvars -except saveFigures
clc

% loadsol = sprintf('/Volumes/MELANOMAII/Data/CriteriaAnalysis/5nodes/sol10005_1_92');
loadsol = sprintf('./Data/CriteriaAnalysis/5nodes/sol10005_1_92');
load(loadsol);
param = 5;

%histogram
pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{param}.samp(1,:),'Facecolor',[135 170 222]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1,'Binwidth',10);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([0, 0.5, 1])
yticklabels({'0', '50', '100'})
xticks([0, 250])
xlabel('gene product')
ylabel('% of cells')
xlim([0,250])
ylim([0 1])

%carpet (from figures before, new)
set(0, 'DefaultFigureRenderer', 'painters');
i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
plot(sol{param}.samp(1,:),ones(length(sol{param}),1),'.','Color',[135 170 222]./255);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca, 'visible', 'off')
set(0, 'DefaultFigureRenderer', 'painters');
xlim([0,250])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure5/Figure5F'])
    print('-dpdf',['./Figures/Figure5/Figure5F'])
end

%% Figure5G

clearvars -except saveFigures
clc

loadsol = sprintf('./Data/CriteriaAnalysis/5nodes/sol10005_10_92');
load(loadsol);
param = 5;

%histogram
pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{param}.samp(1,:),'Facecolor',[135 170 222]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1,'Binwidth',10);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([0, 0.5, 1])
yticklabels({'0', '50', '100'})
xticks([0, 250])
xlabel('gene product')
ylabel('% of cells')
xlim([0,250])
ylim([0 1])

%carpet (from figures before, new)
set(0, 'DefaultFigureRenderer', 'painters');
i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
plot(sol{param}.samp(1,:),ones(length(sol{param}),1),'.','Color',[135 170 222]./255);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca, 'visible', 'off')
set(0, 'DefaultFigureRenderer', 'painters');
xlim([0,250])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure5/Figure5G'])
    print('-dpdf',['./Figures/Figure5/Figure5G'])
end

%% Figure5H

clearvars -except saveFigures
clc

loadPhixerData

thresPhi = dataforplots(:,1);
val = dataforplots(:,2:end);

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
plot(thresPhi,val(:,1:2),'-','Color',[211,95,95]./255,'Linewidth',1)
hold on
plot(thresPhi,val(:,3:end),'-','Color','k','Linewidth',1)
set(gca,'FontSize',12)
yticks([0, 125, 250])
xticks([0.4, 0.7])
xlabel('edge weight threshold')
ylabel('no. edges')
xlim([0.4, 0.7])
ylim([0 250])
box off
set(gca,'linewidth',2)

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure5/Figure5H'])
    print('-dpdf',['./Figures/Figure5/Figure5H'])
end

%% Figure5I

clearvars -except saveFigures
clc

ImportNoDrug2Bootstrapped
ImportDrug5Bootstrapped

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(nodrug2,'BinWidth',5,'Facecolor',[211,95,95]./255,'Normalization','probability',...
    'Edgecolor','none','FaceAlpha', 1);
hold on
histogram(fourweeks5,'BinWidth',5,'Facecolor',[50,50,50]./255,'Normalization','probability',...
    'Edgecolor','none','FaceAlpha', 1);
xlim([0,220])
ylim([0,1])
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([0, 1])
xticks([0, 110, 220])
xlabel('no. edges')
ylabel('frequency')
xlim([0,220])
ylim([0 1])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure5/Figure5I'])
    print('-dpdf',['./Figures/Figure5/Figure5I'])
end
