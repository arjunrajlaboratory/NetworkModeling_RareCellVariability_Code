
addpath(genpath(pwd));

saveFigures = true;

%% SupplementeryFigure4B 

% load('/Volumes/MELANOMA/Data/Data1000.mat')
load('./Data/Data1000.mat')
n_species = 5;
tcell = 1000;    %time per 'cell'
loadsoutpar = sprintf('./Data/Translation/S_outparTransFaster5_3_1');
a = 10;
b = 10;

load(loadsoutpar)
thres = Data1000(968,1)^2/Data1000(968,2)^2*Data1000(968,8)*0.8*a/b;
 
i = 1;
param = 1;
threshold = thres;

[maxjackpot_sol,desc_sol,rightskewed_sol,unimodal_sol,samp_sol,samp_time_sol,rand_time_sol] ...
    = analyzeQualTrans_revision(n_species,tcell,threshold,S_outpar{i});

sol{i}.maxjackpot = maxjackpot_sol;
sol{i}.desc = desc_sol;
sol{i}.rightskewed = rightskewed_sol;
sol{i}.unimodal = unimodal_sol;
sol{i}.samp = samp_sol;
sol{i}.samp_time = samp_time_sol;
sol{i}.time = rand_time_sol;
    
if  sol{i}.maxjackpot == 1
    if sol{i}.desc == 1
        if sol{i}.rightskewed == 1
            if sol{i}.unimodal == 1
                rare_par = i;
            end
        end
    end
end
   

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
%number of highly expressed genes
histogram(sol{1}.samp_time,'Facecolor',[0,0,0]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1);
set(gca,'linewidth',2)
box off
ylim([0,1])
yticks([0, 0.5, 1])
yticklabels({'0', '50', '100'})
xlim([-0.5,5.5])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
xlabel('no. highly expressed genes')
ylabel('% of cells')

%expression distribution
j = 2;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
histogram(sol{1}.samp(1,:),'Facecolor',[135 170 222]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1, 'Binwidth',50);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 1400])
xlabel('gene product')
xlim([0,1400])
ylim([0 1])

%carpet
j = 3;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
plot(sol{1}.samp(1,:),ones(1000,1),'.','Color',[135 170 222]./255);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca, 'visible', 'off')
set(0, 'DefaultFigureRenderer', 'painters');
xlim([0,1400])


if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS4/FigureS4B'])
    print('-dpdf',['./Figures/FigureS4/FigureS4B'])
end

%% SupplementeryFigure4E

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);

load('./Data/Multiplicative/sol1000Mult_5_3_97')
%number of highly expressed genes
histogram(sol{8}.samp_time,'Facecolor',[0,0,0]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1);
set(gca,'linewidth',2)
box off
ylim([0,1])
yticks([0, 0.5, 1])
yticklabels({'0', '50', '100'})
xlim([-0.5,5.5])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
xlabel('no. highly expressed genes')
ylabel('% of cells')

%expression distribution
j = 2;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
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
j = 3;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
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
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS4/FigureS4E'])
    print('-dpdf',['./Figures/FigureS4/FigureS4E'])
end

%% SupplementeryFigure4F

clearvars -except saveFigures
clc

for irep = 1:3
    
    loadCom = sprintf('./Data/Com%d',irep);
    load(loadCom)
    
    count = 1;
    for n_species = [2,3,5,8]
        N(count) = length(find(Com(:,2) == n_species & Com(:,7)>Com(:,8)));
        count = count +1;
    end

    %less stringent
    count = 1;
    N_less1 = find(Com(:,1) == 0);
    for n_species = [2,3,5,8]
        N_less2(count) = length(find(Com(N_less1,2) == n_species));
        count = count +1;
    end
    N_less = N_less2+N;
    
    R(irep,:) = N;
    R_less(irep,:) = N_less;
    
end

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
boxplot([R_less(:,1),R_less(:,2), R_less(:,3), R_less(:,4)],'Colors', [0,0,0]./255)
ylim([0,170])
xlabel('network size')
xticklabels({'2','3','5','8'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
ylabel('simulations')
yticks([0,85,170])
box off
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
set(gca,'linewidth',2)


if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS4/FigureS4F'])
    print('-dpdf',['./Figures/FigureS4/FigureS4F'])
end

%% SupplementeryFigure4G

clearvars -except saveFigures
clc

loadCom = sprintf('./Data/Com%d',1);
load(loadCom)

N_less1 = find(Com(:,1) == 0);
for inet = 1:10
    R(inet) = length(find(Com(:,2) == 5 & Com(:,3) == inet & Com(:,7)>Com(:,8)));
    R(inet) = R(inet)+length(find(Com(N_less1,2) == 5 & Com(N_less1,3) == inet));
end

% load('/Volumes/MELANOMA/Data/M_iso5')
load('./Data/M_iso5')
for jnet = 1:length(M_iso)
    C(jnet) = sum(M_iso{jnet}(:,1));
end

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);

plot(C,R,'.','Markersize',20,'Color', [211,95,95]./255)
ylabel('simulations')
xlabel('connectivity')
xticks([1,2,3,4,5])
yticks([0,20,40])
xlim([0.5,5.5])
ylim([0,40])
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca,'linewidth',2)

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS4/FigureS4G'])
    print('-dpdf',['./Figures/FigureS4/FigureS4G'])
end

%% SupplementeryFigure4H

clearvars -except saveFigures
clc

loadCom = sprintf('./Data/Com%d',1);
load(loadCom)

A = find(Com(:,7)>Com(:,8));
for inum = 1:1000
    ParFreq(inum) = length(find(Com(A,6) == inum));
end

N_less1 = find(Com(:,1) == 0);
A = find(Com(:,7)>Com(:,8));
for inum = 1:1000
    ParFreq(inum) = length(find(Com(A,6) == inum));
    ParFreq(inum) = ParFreq(inum) + length(find(Com(N_less1,1) == 0 & Com(N_less1,6) == inum));
end

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);

histogram(ParFreq(ParFreq<15)/96*100,'FaceColor',[0.2, 0.2, 0.2],'Binwidth',1,...
    'Edgecolor','none')
hold on
histogram(ParFreq(ParFreq>=15)/96*100,'FaceColor',[255, 127, 42]./255','Edgecolor','none','Binwidth',1)
alpha(1)
hold on 
plot([15,15],[0,920],':','Color','k','Linewidth',3)
xlim([0,50])
ylim([0,940])
xlabel('% of simulations')
ylabel('Counts')
set(gca,'linewidth',2)
box off
xticks([0,25,50])
yticks([0, 10, 20, 30, 930, 940])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
breakyaxis([30, 930]);

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS4/FigureS4H'])
    print('-dpdf',['./Figures/FigureS4/FigureS4H'])
end


%% SupplementaryFigure4I

clearvars -except saveFigures
clc

% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Data1000')

R_an = Data1000(:,[1:6,8]);

%from previous analysis (Figure2G) we find 8 rare coordinated high parameter sets
I = [26 92 133 544 702 915 968];

DecisionTree(I,R_an)

%x3 - rate r_on
%x5 - rate r_add
%x6 - rate r_off

%0 - no rare coordinated high parameter set
%1 - rare coordinated high parameter set

%% SupplementeryFigure4J

clearvars -except saveFigures

for irep = 1:3
    
    loadCom = sprintf('./data/Com%d',irep);
    load(loadCom)
    
    count = 1;
    for n_species = [2,3,5,8]
        N(count) = length(find(Com(:,2) == n_species & Com(:,7)>Com(:,8)));
        count = count +1;
    end

    %less stringent
    count = 1;
    N_less1 = find(Com(:,1) == 0);
    for n_species = [2,3,5,8]
        N_less2(count) = length(find(Com(N_less1,2) == n_species));
        count = count +1;
    end
    N_less = N_less2+N;
    
    R(irep,:) = N;
    R_less(irep,:) = N_less;
    
end

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);

boxplot([R(:,1),R(:,2), R(:,3), R(:,4)],'Colors', [0,0,0]./255)
ylim([0,140])
xlabel('network size')
xticklabels({'2','3','5','8'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
ylabel('simulations')
yticks([0,70,140])
box off
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
set(gca,'linewidth',2)

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS4/FigureS4J'])
    print('-dpdf',['./Figures/FigureS4/FigureS4J'])
end

%% SupplementeryFigure4K

clearvars -except saveFigures
clc

loadCom = sprintf('./Data/Com%d',1);
load(loadCom)

for inet = 1:10
    R(inet) = length(find(Com(:,2) == 5 & Com(:,3) == inet & Com(:,7)>Com(:,8)));
end

% load('/Volumes/MELANOMA/Data/M_iso5')
load('./Data/M_iso5')
for jnet = 1:length(M_iso)
    C(jnet) = sum(M_iso{jnet}(:,1));
end

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);

plot(C,R,'.','Markersize',20,'Color', [211,95,95]./255)
ylabel('simulations')
xlabel('connectivity')
xticks([1,2,3,4,5])
yticks([0,15,30])
xlim([0.5,5.5])
ylim([0,30])
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca,'linewidth',2)

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS4/FigureS4K'])
    print('-dpdf',['./Figures/FigureS4/FigureS4K'])
end

%% SupplementeryFigure4L

clearvars -except saveFigures
clc

loadCom = sprintf('./Data/Com%d',1);
load(loadCom)

A = find(Com(:,7)>Com(:,8));
for inum = 1:1000
    ParFreq(inum) = length(find(Com(A,6) == inum));
end

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);

histogram(ParFreq(ParFreq<15)/96*100,'FaceColor',[0.2, 0.2, 0.2],'Binwidth',1,...
    'Edgecolor','none')
hold on
histogram(ParFreq(ParFreq>=15)/96*100,'FaceColor',[255, 127, 42]./255','Edgecolor','none','Binwidth',1)
alpha(1)
hold on 
plot([15,15],[0,970],':','Color','k','Linewidth',2)
xlim([0,50])
ylim([0,970])
xlabel('% of simulations')
ylabel('Counts')
set(gca,'linewidth',2)
box off
xticks([0,25,50])
yticks([0, 10, 20, 30, 960, 970])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
breakyaxis([30, 960]);

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS4/FigureS4L'])
    print('-dpdf',['./Figures/FigureS4/FigureS4L'])
end

%% SupplementaryFigure4M

clearvars -except saveFigures
clc

% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Data1000')

R_an = Data1000(:,[1:6,8]);

%from previous analysis (Figure2G) we find 8 rare coordinated high parameter sets
I = [26 92 133 544 915 968];

DecisionTree(I,R_an)

%x3 - rate r_on
%x5 - rate r_add
%x6 - rate r_off

%0 - no rare coordinated high parameter set
%1 - rare coordinated high parameter set

