
addpath(genpath(pwd));
% addpath(genpath('/Volumes/MELANOMAII'));

saveFigures = true;

%% SupplementeryFigure7A

clearvars -except saveFigures
clc

% load('/Volumes/MELANOMA/Data/BurstAnalysis/3nodes/EnterTimes/sol_EnterJack3_2_1')
load('./Data/BurstAnalysis/3nodes/EnterTimes/sol_EnterJack3_2_1')

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);

h = histfit(sol.Data{1}(:,1),10,'exponential');
yt = get(gca, 'YTick');
set(h(1),'facecolor',[100,100,100]./255, 'edgecolor', 'none', 'facealpha', 1);
set(h(2), 'color','k');
xlim([0, 300])
ylim([0, (length(sol.Data{1}(:,1))-sum(isnan(sol.Data{1}(:,1))))*0.65])
xticks([0, 300])
yticks([0, (length(sol.Data{1}(:,1))-sum(isnan(sol.Data{1}(:,1))))*0.65])
box off
xlabel('t_{ent}'),
ylabel('probability')
yticklabels({'0','0.65'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS7/FigureS7A'])
    print('-dpdf',['./Figures/FigureS7/FigureS7A'])
end

%% SupplementeryFigure7B

clearvars -except saveFigures
clc

% load('/Volumes/MELANOMA/Data/BurstAnalysis/3nodes/ExitTimes/sol_ExitJack3_2_2')
load('./Data/BurstAnalysis/3nodes/ExitTimes/sol_ExitJack3_2_2')
pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
h = histfit(sol.Data{5}(:,1),10,'exponential');
y = get(gca, 'YTick');
set(h(1),'facecolor',[100,100,100]./255, 'edgecolor', 'none', 'facealpha', 1);
set(h(2), 'color','k');
xlim([0, 210])
ylim([0, (length(sol.Data{5}(:,1))-sum(isnan(sol.Data{5}(:,1))))])
xticks([0, 210])
yticks([0, (length(sol.Data{5}(:,1))-sum(isnan(sol.Data{5}(:,1))))])
box off
xlabel('t_{exit}'),
ylabel('probability')
yticklabels({'0','1'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS7/FigureS7B'])
    print('-dpdf',['./Figures/FigureS7/FigureS7B'])
end

