
addpath(genpath(pwd));

saveFigures = true;

%% Figure4A

clearvars -except saveFigures
clc

load('./Data/Simulations/3nodes/S_outpar3_2_1')
load('./Data/Data100')
thres = Data100(:,1)./Data100(:,2).*Data100(:,5)*0.8;
param = 88;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.5;
ylen = 0.25;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(856000:856500,S_outpar{param}(1:3,856000:856500),'LineWidth',2);
hold on 
plot(856000:856500,repmat(thres(param),1,length(856000:856500)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'}, {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255});
xlim([856000, 856500])
xticks([856000 856500])
xticklabels({'0', '500'})
yticks([0, 300])
ylim([0,300])
yticklabels({'0', '3'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
ylabel('gene product (10^2)')
xlabel('time')
box off

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure4/Figure4A'])
    print('-dpdf',['./Figures/Figure4/Figure4A'])
end

%% Figure4B

%BurstFraction

% load('/Volumes/MELANOMAII/Data/BurstAnalysis/bursts_high')
% load('/Volumes/MELANOMAII/Data/BurstAnalysis/length_high')
% load('/Volumes/MELANOMAII/Data/BurstAnalysis/bursts_base')
% load('/Volumes/MELANOMAII/Data/BurstAnalysis/length_base')

load('./Data/BurstAnalysis/bursts_high')
load('./Data/BurstAnalysis/length_high')
load('./Data/BurstAnalysis/bursts_base')
load('./Data/BurstAnalysis/length_base')

[h,p,ks2stat] = kstest2((bursts_high./length_high),(bursts_base./length_base));
p

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.4;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
[h,L,MX,MED]=violin([(bursts_high./length_high)',(bursts_base./length_base)'], 'facecolor', [100,100,100;100,100,100]./255);
xticks([1,2])
yticks([0,1])
xticklabels({'high time-region','baseline time-region'})
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
box off
set(gca,'linewidth',2)
ylabel('burst fraction')

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure4/Figure4B'])
    print('-dpdf',['./Figures/Figure4/Figure4B'])
end

%% Figure4C

%BurstFreq

% load('/Volumes/MELANOMAII/Data/BurstAnalysis/num_high')
% load('/Volumes/MELANOMAII/Data/BurstAnalysis/lengthnum_high')
% load('/Volumes/MELANOMAII/Data/BurstAnalysis/num_base')
% load('/Volumes/MELANOMAII/Data/BurstAnalysis/lengthnum_base')

load('./Data/BurstAnalysis/num_high')
load('./Data/BurstAnalysis/lengthnum_high')
load('./Data/BurstAnalysis/num_base')
load('./Data/BurstAnalysis/lengthnum_base')

[h,p,ks2stat] = kstest2((num_high./lengthnum_high),(num_base./lengthnum_base));
p

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.4;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
[h,L,MX,MED]=violin([(num_high./lengthnum_high)',(num_base./lengthnum_base)'], 'facecolor', [100,100,100;100,100,100]./255);
xticks([1,2])
yticks([0,0.1])
xticklabels({'high time-region','baseline time-region'})
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
box off
set(gca,'linewidth',2)
ylabel('burst frequency')
ylim([-0.03,0.1])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure4/Figure4C'])
    print('-dpdf',['./Figures/Figure4/Figure4C'])
end

%% Figure 4D 

%BurstLengthsInd

% load('/Volumes/MELANOMAII/Data/BurstAnalysis/length_L_32_dep')
% load('/Volumes/MELANOMAII/Data/BurstAnalysis/length_L_32')
load('./Data/BurstAnalysis/length_L_32_dep')
load('./Data/BurstAnalysis/length_L_32')

for i = 1:26
    Sdep(i) = sum(length_L_32_dep{i});
    S(i) = sum(length_L_32{i});
    Ndep(i) = length(length_L_32_dep{i});
    N(i) = length(length_L_32{i});
end


pos = [0.1,0.52,0.54,0.76];
xlen = 0.4;
ylen = 0.4;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
[h,L,MX,MED]=violin((Ndep./N)', 'facecolor', [100,100,100;100,100,100]./255);
hold on
plot(ones(length((Ndep./N)),1)+1/50*randn(length((Ndep./N)),1),(Ndep./N),'.','Markersize', 20, 'Color', 'k')
hold on 
plot([0,2],[1,1],':','Color','k','Linewidth',2)
xticks([])
% yticks([0,100])
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
box off
set(gca,'linewidth',2)
xlabel('number of high states')
ylabel('fold change')
ylim([-2,12])


i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
[h,L,MX,MED]=violin((Sdep./S)', 'facecolor', [100,100,100;100,100,100]./255);
hold on
plot(ones(length((Sdep./S)),1)+1/50*randn(length((Sdep./S)),1),(Sdep./S),'.','Markersize', 20, 'Color', 'k')
hold on 
plot([0,2],[1,1],':','Color','k','Linewidth',2)
xticks([])
% yticks([0,100])
xticklabels({'network','indepedent'})
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
box off
set(gca,'linewidth',2)
xlabel('total time in high states')
ylim([-2,12])
yticks([])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure4/Figure4D'])
    print('-dpdf',['./Figures/Figure4/Figure4D'])
end

%% Figure4E

clearvars -except saveFigures
clc

% load('/Volumes/MELANOMA/Data/BurstAnalysis/3nodes/EnterBurstDuration/sol_BurstsOnInit3_2_2')
load('./Data/BurstAnalysis/3nodes/EnterBurstDuration/sol_BurstsOnInit3_2_2')

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram([M.BurstOnJInit_Data{4,1},M.BurstOnJInit_Data{4,2},M.BurstOnJInit_Data{4,3}],...
    'Binwidth',10,'Normalization','probability','EdgeColor','none','FaceColor', [100,100,100]./255,...
    'FaceAlpha',1)
hold on
histogram([M.BurstOnJNInit_Data{4,1},M.BurstOnJNInit_Data{4,2},M.BurstOnJNInit_Data{4,3}],...
    'Binwidth',10,'Normalization','probability','EdgeColor','none','FaceColor',[0,0,0],...
    'FaceAlpha',1)
xlim([0,260])
ylim([0,0.7])
box off
xlabel('duration'),
ylabel('probability')
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
yticks([0,0.7])
xticks([0, 260])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure4/Figure4E'])
    print('-dpdf',['./Figures/Figure4/Figure4E'])
end

%% Figure4F

clearvars -except saveFigures
clc

% load('/Volumes/MELANOMA/Data/BurstAnalysis/3nodes/ExitBurstDuration/sol_BurstsTerm3_2_2')
load('./Data/BurstAnalysis/3nodes/ExitBurstDuration/sol_BurstsTerm3_2_2')

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram([M.BurstOnJExit_Data{4,1},M.BurstOnJExit_Data{4,2},M.BurstOnJExit_Data{4,3}],...
    'Binwidth',10,'Normalization','probability','EdgeColor','none','FaceColor', [150,150,150]./255,...
    'FaceAlpha',1)
hold on
histogram([M.BurstOnJNExit_Data{4,1},M.BurstOnJNExit_Data{4,2},M.BurstOnJNExit_Data{4,3}],...
    'Binwidth',10,'Normalization','probability','EdgeColor','none','FaceColor',[204, 204, 204]./255,...
    'FaceAlpha',1)
xlim([0,260])
ylim([0,0.7])
box off
xlabel('duration'),
ylabel('probability')
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
yticks([0,0.7])
xticks([0, 260])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure4/Figure4F'])
    print('-dpdf',['./Figures/Figure4/Figure4F'])
end

%% Figure4G

clearvars -except saveFigures
clc

I = [];
T = [];
IEval = [];
TEval = [];


for n_species = [2,3,5,8]
%     loadInit = sprintf('/Volumes/MELANOMAII/Data/EvalRatioInit%d',n_species);
    loadInit = sprintf('./Data/EvalRatioInit%d',n_species);
    load(loadInit)
%     loadEvalInit = sprintf('/Volumes/MELANOMAII/Data/EvalInit%d',n_species);
    loadEvalInit = sprintf('./Data/EvalInit%d',n_species);
    load(loadEvalInit)
%     loadTerm = sprintf('/Volumes/MELANOMAII/Data/EvalRatioTerm%d',n_species);
    loadTerm = sprintf('./Data/EvalRatioTerm%d',n_species);
    load(loadTerm)
%     loadEvalTerm = sprintf('/Volumes/MELANOMAII/Data/EvalTerm%d',n_species);
    loadEvalTerm = sprintf('./Data/EvalTerm%d',n_species);
    load(loadEvalTerm)
    
    I = [I, ROnInit];
    T = [T, ROnTerm];
    
    IEval = [IEval;MAll];
    TEval = [TEval;MAllOn];
    
end

I0 = I(IEval==0);
T1 = (T(TEval==1));
T0 = (T(TEval==0));

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.4;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
[h,L,MX,MED]=violin([I(1:488)',T(~isnan(T(1:491)))'], 'facecolor', [150,150,150;150,150,150]./255)
xticks([1,2])
yticks([0,30,60])
xticklabels({'entry','exit'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
box off
set(gca,'linewidth',2)
ylabel('mean burst duration ratio')

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure4/Figure4G'])
    print('-dpdf',['./Figures/Figure4/Figure4G'])
end

%% Figure4H

clearvars -except saveFigures
clc

% load('/Volumes/MELANOMA/Data/BurstAnalysis/3nodes/EnterTimes/sol_EnterJack3_2_2')
load('./Data/BurstAnalysis/3nodes/EnterTimes/sol_EnterJack3_2_2')

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = histfit(sol.Data{4}(:,1),10,'exponential');
yt = get(gca, 'YTick');
set(h(1),'facecolor',[100,100,100]./255, 'edgecolor', 'none', 'facealpha', 1);
set(h(2), 'color','k');
xlim([0, 150])
ylim([0, (length(sol.Data{4}(:,1))-sum(isnan(sol.Data{4}(:,1))))*0.65])
xticks([0, 150])
yticks([0, (length(sol.Data{4}(:,1))-sum(isnan(sol.Data{4}(:,1))))*0.65])
box off
xlabel('t_{ent}'),
ylabel('probability')
yticklabels({'0','0.65'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure4/Figure4H'])
    print('-dpdf',['./Figures/Figure4/Figure4H'])
end

%% Figure4I

clearvars -except saveFigures
clc

% load('/Volumes/MELANOMA/Data/BurstAnalysis/3nodes/ExitTimes/sol_ExitJack3_2_2')
load('./Data/BurstAnalysis/3nodes/ExitTimes/sol_ExitJack3_2_2')

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = histfit(sol.Data{1}(:,1),10,'exponential');
set(h(1),'facecolor',[100,100,100]./255, 'edgecolor', 'none', 'facealpha', 1);
set(h(2), 'color','k');
xlim([0, 150])
ylim([0, (length(sol.Data{1}(:,1))-sum(isnan(sol.Data{1}(:,1))))*0.65])
xticks([0, 150])
yticks([0, (length(sol.Data{1}(:,1))-sum(isnan(sol.Data{1}(:,1))))*0.65])
box off
xlabel('t_{exit}'),
ylabel('probability')
yticklabels({'0','0.65'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure4/Figure4I'])
    print('-dpdf',['./Figures/Figure4/Figure4I'])
end
