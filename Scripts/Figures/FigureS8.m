

addpath(genpath(pwd));
% addpath(genpath('/Volumes/MELANOMAII'));

saveFigures = true;

%% SupplementaryFigure8A

clearvars -except saveFigures
clc

ImportNoDrug1Bootstrapped;
ImportDrug1Bootstrapped;
ImportDrug2Bootstrapped;
ImportDrug3Bootstrapped;
ImportDrug4Bootstrapped;
ImportDrug6Bootstrapped;
ImportDrug7Bootstrapped;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);

histogram(nodrug1,'BinWidth',5,'Facecolor',[211,95,95]./255,'Normalization','probability',...
    'Edgecolor','none','FaceAlpha', 0.8);
hold on
histogram(fourweeks1,'BinWidth',5,'Facecolor',[50,50,50]./255,'Normalization','probability',...
    'Edgecolor','none','FaceAlpha', 0.8);
hold on
histogram(fourweeks2,'BinWidth',5,'Facecolor',[50,50,50]./255,'Normalization','probability',...
    'Edgecolor','none','FaceAlpha', 0.8);
hold on
histogram(fourweeks3,'BinWidth',5,'Facecolor',[50,50,50]./255,'Normalization','probability',...
    'Edgecolor','none','FaceAlpha', 0.8);
hold on
histogram(fourweeks4,'BinWidth',5,'Facecolor',[50,50,50]./255,'Normalization','probability',...
    'Edgecolor','none','FaceAlpha', 0.8);
hold on
histogram(fourweeks6,'BinWidth',5,'Facecolor',[50,50,50]./255,'Normalization','probability',...
    'Edgecolor','none','FaceAlpha', 0.8);
hold on
histogram(fourweeks7,'BinWidth',5,'Facecolor',[50,50,50]./255,'Normalization','probability',...
    'Edgecolor','none','FaceAlpha', 0.8);
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
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS8/FigureS8A'])
    print('-dpdf',['./Figures/FigureS8/FigureS8A'])
end


%% SupplementaryFigure8B

clearvars -except saveFigures
clc

ImportNoDrug1Control
ImportNoDrug2Control
ImportDrug1Control
ImportDrug2Control
ImportDrug3Control
ImportDrug4Control
ImportDrug5Control
ImportDrug6Control
ImportDrug7Control

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
histogram(controlnodrug1,'BinWidth',0.05,'Facecolor',[211,95,95]./255,'Normalization','probability',...
    'Edgecolor','none','FaceAlpha', 1);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([0, 0.5])
xticks([0, 0.6])
xlabel('edge weight')
ylabel('frequency')
xlim([0,0.6])
ylim([0 0.5])

j = 2;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
histogram(controlnodrug2,'BinWidth',0.05,'Facecolor',[211,95,95]./255,'Normalization','probability',...
    'Edgecolor','none','FaceAlpha', 1);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 0.6])
xlabel('edge weight')
xlim([0,0.6])
ylim([0 0.5])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS8/FigureS8Bnodrug'])
    print('-dpdf',['./Figures/FigureS8/FigureS8Bnodrug'])
end

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
histogram(controlfourweeks1cluster1,'BinWidth',0.05,'Facecolor',[50,50,50]./255,'Normalization','probability',...
    'Edgecolor','none','FaceAlpha', 1);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([0, 0.5])
xticks([0, 0.6])
xlabel('edge weight')
ylabel('frequency')
xlim([0,0.6])
ylim([0 0.5])

j = 2;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
histogram(controlfourweeks1cluster2,'BinWidth',0.05,'Facecolor',[50,50,50]./255,'Normalization','probability',...
    'Edgecolor','none','FaceAlpha', 1);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 0.6])
xlabel('edge weight')
xlim([0,0.6])
ylim([0 0.5])

j = 3;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
histogram(controlfourweeks1cluster3,'BinWidth',0.05,'Facecolor',[50,50,50]./255,'Normalization','probability',...
    'Edgecolor','none','FaceAlpha', 1);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 0.6])
xlabel('edge weight')
xlim([0,0.6])
ylim([0 0.5])

j = 4;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
histogram(controlfourweeks1cluster4,'BinWidth',0.05,'Facecolor',[50,50,50]./255,'Normalization','probability',...
    'Edgecolor','none','FaceAlpha', 1);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 0.6])
xlabel('edge weight')
xlim([0,0.6])
ylim([0 0.5])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS8/FigureS8BdrugI'])
    print('-dpdf',['./Figures/FigureS8/FigureS8BdrugI'])
end

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
histogram(controlfourweeks2cluster1,'BinWidth',0.05,'Facecolor',[50,50,50]./255,'Normalization','probability',...
    'Edgecolor','none','FaceAlpha', 1);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([0, 0.5])
xticks([0, 0.6])
xlabel('edge weight')
ylabel('frequency')
xlim([0,0.6])
ylim([0 0.5])

j = 2;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
histogram(controlfourweeks2cluster2,'BinWidth',0.05,'Facecolor',[50,50,50]./255,'Normalization','probability',...
    'Edgecolor','none','FaceAlpha', 1);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 0.6])
xlabel('edge weight')
xlim([0,0.6])
ylim([0 0.5])

j = 3;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
histogram(controlfourweeks2cluster3,'BinWidth',0.05,'Facecolor',[50,50,50]./255,'Normalization','probability',...
    'Edgecolor','none','FaceAlpha', 1);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 0.6])
xlabel('edge weight')
xlim([0,0.6])
ylim([0 0.5])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS8/FigureS8BdrugII'])
    print('-dpdf',['./Figures/FigureS8/FigureS8BdrugII'])
end

%% SupplementaryFigure8C

%original networks can be found in /Volumes/MELANOMAII/Data/Networks

