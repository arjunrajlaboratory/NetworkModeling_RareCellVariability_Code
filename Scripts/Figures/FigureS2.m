
addpath(genpath(pwd));

saveFigures = true;

%% SupplementaryFigure2A

%0 - stabel low
%1 - uncorrelated transient high
%2 - rare transient correlaetd high
%3 - stably high

class2 = load('./Data/classAll2');
class3 = load('./Data/classAll3');
class5 = load('./Data/classAll5');
class8 = load('./Data/classAll8');

C = [class2.class(:,4);class3.class(:,4);class5.class(:,4);class8.class(:,4)];

Num0 = sum(C == 0)/length(C)
Num1 = sum(C == 1)/length(C)
Num2 = sum(C == 2)/length(C)
Num3 = sum(C == 3)/length(C)

%% SupplementaryFigure2B

clearvars -except saveFigures
clc

%SupplementaryFigure2B top row

load('./Data/CriteriaAnalysis/2nodes/sol10002_1_6')
param = 44;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{param}.samp_time,'Facecolor',[0,0,0]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1);
set(gca,'linewidth',2)
box off
ylim([0,1])
yticks([0, 0.5, 1])
yticklabels({'0', '50', '100'})
xlim([-0.5,2.5])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
xlabel('no. highly expressed genes')
ylabel('% of cells')

i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{param}.samp(1,:),'Facecolor',[135 170 222]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1, 'Binwidth',50)
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 750])
xlabel('gene product')
xlim([0,750])
ylim([0 1])

i = 3;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
plot(sol{param}.samp(1,:),ones(1000,1),'.','Color',[135 170 222]./255);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca, 'visible', 'off')
set(0, 'DefaultFigureRenderer', 'painters');
xlim([0,750])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS2/FigureS2Btop'])
    print('-dpdf',['./Figures/FigureS2/FigureS2Btop'])
end

%SupplementaryFigure2B middle row

clearvars -except saveFigures
clc

load('./Data/CriteriaAnalysis/3nodes/sol10003_2_6')
param = 44;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{param}.samp_time,'Facecolor',[0,0,0]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1);
set(gca,'linewidth',2)
box off
ylim([0,1])
yticks([0, 0.5, 1])
yticklabels({'0', '50', '100'})
xlim([-0.5,3.5])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
xlabel('no. highly expressed genes')
ylabel('% of cells')

i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{param}.samp(1,:),'Facecolor',[135 170 222]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1, 'Binwidth',50)
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 750])
xlabel('gene product')
xlim([0,750])
ylim([0 1])

i = 3;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
plot(sol{param}.samp(1,:),ones(1000,1),'.','Color',[135 170 222]./255);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca, 'visible', 'off')
set(0, 'DefaultFigureRenderer', 'painters');
xlim([0,750])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS2/FigureS2Bmiddle'])
    print('-dpdf',['./Figures/FigureS2/FigureS2Bmiddle'])
end

%SupplementaryFigure2B bottom row

clearvars -except saveFigures
clc

load('./Data/CriteriaAnalysis/5nodes/sol10005_1_55')
param = 4;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{param}.samp_time,'Facecolor',[0,0,0]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1);
set(gca,'linewidth',2)
box off
ylim([0,1])
yticks([0, 0.5, 1])
yticklabels({'0', '50', '100'})
xticks(0:5)
xlim([-0.5,5.5])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
xlabel('no. highly expressed genes')
ylabel('% of cells')

i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{param}.samp(1,:),'Facecolor',[135 170 222]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1, 'Binwidth',50)
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 750])
xlabel('gene product')
xlim([0,750])
ylim([0 1])

i = 3;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
plot(sol{param}.samp(1,:),ones(1000,1),'.','Color',[135 170 222]./255);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca, 'visible', 'off')
set(0, 'DefaultFigureRenderer', 'painters');
xlim([0,750])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS2/FigureS2Bbottom'])
    print('-dpdf',['./Figures/FigureS2/FigureS2Bbottom'])
end

%% SupplementaryFigure2C

clearvars -except saveFigures
clc

load('./Data/CriteriaAnalysis/5nodes/sol10005_1_55')
param = 4;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{param}.samp(2,:),'Facecolor',[22 45 80]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1, 'Binwidth',50)
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([0, 0.5, 1])
yticklabels({'0', '50', '100'})
xticks([])
ylabel('% of cells')
xlim([0,750])
ylim([0 1])
xticks([0, 750])
xlabel('gene product')

i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{param}.samp(3,:),'Facecolor',[95 141 211]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1, 'Binwidth',50)
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 750])
xlabel('gene product')
xlim([0,750])
ylim([0 1])

i = 3;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{param}.samp(4,:),'Facecolor',[44 90 160]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1, 'Binwidth',50)
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 750])
xlabel('gene product')
xlim([0,750])
ylim([0 1])

i = 4;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{param}.samp(5,:),'Facecolor',[33 68 120]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1, 'Binwidth',50)
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 750])
xlabel('gene product')
xlim([0,750])
ylim([0 1])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS2/FigureS2C'])
    print('-dpdf',['./Figures/FigureS2/FigureS2C'])
end

%% SupplementaryFigure2E

%10nodes

loadsoutpar = sprintf('./Data/10nodes/S_outpar1210_1_1');
load(loadsoutpar)

loadsol = sprintf('./Data/10nodes/sol1000_10_1');
load(loadsol)
    
% load('/Volumes/MELANOMA/Data/Data1000.mat')
load('./Data/Data1000.mat')
PAll = [26 92 133 183 544 702 915 968 504 889 125 944];

Data12 = Data1000(PAll,:);
thres = Data12(:,1)./Data12(:,2).*Data12(:,8)*0.8;

%colors?
pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.4;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(0:300,S_outpar{8}(1:10,447100:447400),'LineWidth',2);
set(h, {'color'}, {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255; [44 90 160]./255;...
    [33 68 120]./255; [85 153 255]./255; [135 170 222]./255; [44 90 160]./255; [170 204 255]./255; [0 102 255]./255});
hold on
plot([0,300],[thres(8),thres(8)],':','LineWidth',...
    2, 'Color', 'k')
yticks([0, 350])
yticklabels({'0', '3.5'})
ylim([0,350])
xticks([0,300])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
ylabel('gene product (10^2)')
xlabel('time')
box off

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS2/FigureS2E'])
    print('-dpdf',['./Figures/FigureS2/FigureS2E'])
end

%% SupplementaryFigure2F

loadsol = sprintf('./Data/10nodes/sol1000_10_1');
load(loadsol)

%right panel 'Simulation'
pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
param = 8;
histogram(sol{8}.samp_time,'Facecolor',[0,0,0]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1);
set(gca,'linewidth',2)
box off
ylim([0,1])
yticks([0, 0.5, 1])
yticklabels({'0', '50', '100'})
xlim([-0.5,10.5])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
xlabel('no. highly expressed genes')
ylabel('% of cells')

%right panel 'Simulation'
i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{8}.samp(1,:),'Facecolor',[135 170 222]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1, 'Binwidth',50);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 350])
xlabel('gene product')
xlim([0,350])
ylim([0 1])

%carpet
i = 3;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
plot(sol{8}.samp(1,:),ones(1000,1),'.','Color',[135 170 222]./255);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca, 'visible', 'off')
set(0, 'DefaultFigureRenderer', 'painters');
xlim([0,350])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS2/FigureS2F'])
    print('-dpdf',['./Figures/FigureS2/FigureS2F'])
end

%% SupplementeryFigure2H 

clearvars -except saveFigures
clc

load('./Data/Simulations/AsymmetricArchitecture/S_outpar5_1_1')
% load('/Volumes/MELANOMA/Data/Data100')
load('./Data/Data100')
thres = Data100(:,1)./Data100(:,2).*Data100(:,5)*0.8;
param = 88;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.4;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(539500:540500,S_outpar{param}(1:5,539500:540500),'LineWidth',2);
hold on 
plot(539500:540500,repmat(thres(param),1,length(539500:540500)),':','LineWidth',...
    2, 'Color', 'k')
set(h, {'color'}, {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255;...
    [44 90 160]./255; [33 68 120]./255});
xlim([539500, 540500])
xticks([539500, 540500])
xticklabels({'0', '1000'})
yticks([0, 350])
yticklabels({'0', '3.5'})
ylim([0,350])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
ylabel('gene product (10^2)')
xlabel('time')
box off

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS2/FigureS2H'])
    print('-dpdf',['./Figures/FigureS2/FigureS2H'])
end

%% SupplementeryFigure2I

clearvars -except saveFigures
clc

load('./Data/CriteriaAnalysis/AsymmetricArchitecture/solAsymmetricArchitectures5_1_1')

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{88}.samp_time,'Facecolor',[0,0,0]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1);
set(gca,'linewidth',2)
box off
ylim([0,1])
yticks([0, 0.5, 1])
yticklabels({'0', '50', '100'})
xlim([-0.5,5.5])
xticks([0:5])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
xlabel('no. highly expressed genes')
ylabel('% of cells')

%right panel - histogram 
i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{88}.samp(1,:),'BinWidth',50,...
    'Facecolor',[135 170 222]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1)
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 400])
xlabel('gene product')
xlim([0,400])
ylim([0 1])

%right panel - carpet
i = 3;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
plot(sol{88}.samp(1,:),ones(1000,1),'.','Color',[135 170 222]./255);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca, 'visible', 'off')
set(0, 'DefaultFigureRenderer', 'painters');
xlim([0,400])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS2/FigureS2I'])
    print('-dpdf',['./Figures/FigureS2/FigureS2I'])
end

%% SupplementeryFigure2K

clearvars -except saveFigures
clc

load('./Data/Simulations/AsymmetricParameterSet/S_outpar5_3_1')
load('./Data/Data100Asym')
thres = Data100Asym(:,1:5)./Data100Asym(:,6:10).*Data100Asym(:,21:25)*0.8;
param = 35;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.4;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(290200:291200,S_outpar{param}(1:5,290200:291200),'LineWidth',2);
hold on 
plot(290200:291200,repmat(thres(param,1),1,length(290200:291200)),':','LineWidth',...
    2, 'Color', [135 170 222]./255)
plot(290200:291200,repmat(thres(param,2),1,length(290200:291200)),':','LineWidth',...
    2, 'Color', [22 45 80]./255)
plot(290200:291200,repmat(thres(param,3),1,length(290200:291200)),':','LineWidth',...
    2, 'Color', [95 141 211]./255)
plot(290200:291200,repmat(thres(param,4),1,length(290200:291200)),':','LineWidth',...
    2, 'Color', [44 90 160]./255)
plot(290200:291200,repmat(thres(param,5),1,length(290200:291200)),':','LineWidth',...
    2, 'Color', [33 68 120]./255)
set(h, {'color'}, {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255;...
    [44 90 160]./255; [33 68 120]./255});
xlim([290200, 291200])
xticks([290200 291200])
xticklabels({'0', '1000'})
yticks([0, 700 1400])
yticklabels({'0', '0.7', '1.4'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
ylabel('gene product (10^3)')
xlabel('time')
box off

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS2/FigureS2K'])
    print('-dpdf',['./Figures/FigureS2/FigureS2K'])
end

%% SupplementeryFigure2L

clearvars -except saveFigures
clc

load('./Data/CriteriaAnalysis/AsymmetricParameterSet/sol5_3_1')

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);

histogram(sol{35}.samp_time,'Facecolor',[0,0,0]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1);
set(gca,'linewidth',2)
box off
ylim([0,1])
yticks([0, 0.5, 1])
yticklabels({'0', '50', '100'})
xlim([-0.5,5.5])
xticks([0:5])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
xlabel('no. highly expressed genes')
ylabel('% of cells')

%right panel - histogram 
i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{35}.samp(1,:),'BinWidth',50,...
    'Facecolor',[135 170 222]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1)
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 1000])
xlabel('gene product')
xlim([0,1000])
ylim([0 1])

%right panel - carpet
i = 3;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
plot(sol{35}.samp(1,:),ones(1000,1),'.','Color',[135 170 222]./255);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca, 'visible', 'off')
set(0, 'DefaultFigureRenderer', 'painters');
xlim([0,1000])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS2/FigureS2L'])
    print('-dpdf',['./Figures/FigureS2/FigureS2L'])
end

%% SupplementeryFigure2M

clearvars -except saveFigures
clc

load('./Data/CriteriaAnalysis/AsymmetricParameterSet/sol5_3_1')
pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{35}.samp(2,:),'BinWidth',max(sol{35}.samp(2,:)/20),...
    'Facecolor',[22 45 80]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1)
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([0, 0.5, 1])
yticklabels({'0', '50', '100'})
xticks([0, 1500])
xlabel('gene product')
ylabel('% of cells')
xlim([0,1500])
ylim([0 1])

i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{35}.samp(3,:),'BinWidth',max(sol{35}.samp(3,:)/20),...
    'Facecolor',[95 141 211]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1)
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 17])
xlabel('gene product')
xlim([0,17])
ylim([0 1])

i = 3;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{35}.samp(4,:),'BinWidth',max(sol{35}.samp(4,:)/20),...
    'Facecolor',[44 90 160]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1)
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 200])
xlabel('gene product')
xlim([0,200])
ylim([0 1])

i = 4;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{35}.samp(5,:),'BinWidth',max(sol{35}.samp(5,:)/20),...
    'Facecolor',[33 68 120]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1)
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 850])
xlabel('gene product')
xlim([0,850])
ylim([0 1])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS2/FigureS2M'])
    print('-dpdf',['./Figures/FigureS2/FigureS2M'])
end
