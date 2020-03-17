addpath(genpath(pwd));

saveFigures = true;

%% Figure3A

clearvars -except saveFigures
clc

JP = zeros(1,0);

for n_species = [2,3,5,8]
    
%     loadMiso = sprintf('/Volumes/MELANOMA/Data/M_iso%d',n_species);
    loadMiso = sprintf('./Data/M_iso%d',n_species);
    load(loadMiso);
    
    for inet = 1:length(M_iso)
%         loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/rare_par1000_%d_%d',n_species,n_species,inet);
        loadrare = sprintf('./Data/RareParameters/%dnodes/rare_par1000_%d_%d',n_species,n_species,inet);
        load(loadrare)
        if isempty(rare_par) == 0
            JP = [JP,rare_par];
        end
    end
end

for inum = 1:1000
    ParFreq(inum) = length(find(JP == inum));
end

pos = [0.1,0.32,0.54,0.76];
xlen = 0.6;
ylen = 0.25;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(ParFreq(ParFreq<20)/96*100,'FaceColor',[0.2, 0.2, 0.2],'Binwidth',1,...
    'Edgecolor','none')
hold on
histogram(ParFreq(ParFreq>=20)/96*100,'FaceColor',[255, 127, 42]./255','Edgecolor','none','Binwidth',1)
alpha(1)
hold on 
plot([20,20],[0,920],':','Color','k','Linewidth',3)
xlim([0,50])
ylim([0,920])
xlabel('% of simulations with rare coordinated high states per parameter set')
ylabel('Counts')
set(gca,'linewidth',2)
box off
xticks([0,25,50])
yticks([0, 10, 20, 30, 910, 920])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
breakyaxis([30, 910]);

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure3/Figure3A'])
    print('-dpdf',['./Figures/Figure3/Figure3A'])
end

%% Figure3B

clearvars -except saveFigures
clc

% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Data1000')

R_an = Data1000(:,[1:6,8]);

%from previous analysis (Figure2G) we find 8 rare coordinated high parameter sets
I = [26 92 133 183 544 702 915 968];
doplot = 'yes';

DecisionTree(I,R_an)

%x3 - rate r_on
%x5 - rate r_add
%x6 - rate r_off

%0 - no rare coordinated high parameter set
%1 - rare coordinated high parameter set


%% Figure3C

clearvars -except saveFigures
clc

% load('/Volumes/MELANOMA/Data/Simulations/3nodes/S_outpar1000_3_2_10')
% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Simulations/3nodes/S_outpar1000_3_2_10')
load('./Data/Data1000')

set(0, 'DefaultFigureRenderer', 'painters');
I = [26 92 133 183 544 702 915 968];

pos = [0.1,0.32,0.54,0.76];
xlen = 0.25;
ylen = 0.25;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
plot3(Data1000(:,3),Data1000(:,5),Data1000(:,6),'.','MarkerSize',2,'Color',[100,100,100]./255)
hold on
plot3(Data1000(I,3),Data1000(I,5),Data1000(I,6),'.','MarkerSize',10,'Color',[255,127,42]./255)
hold on
% for i = 1:length(I)
%     p1 = plot3([Data1000(I(i),3),Data1000(I(i),3)],[Data1000(I(i),5),Data1000(I(i),5)],...
%         [zeros(1,1),Data1000(I(i),6)],'--','Color',[255,127,42]./255,'LineWidth',2);
%     p1.Color(4) = 0.25;
%     hold on
%     plot3(Data1000(I(i),3),Data1000(I(i),5),0,'.','MarkerSize',20,'Color',[255,127,42]./255)
%     hold on
% end
xlim([0,0.1])
ylim([0,1])
zlim([0,0.1])
xlabel('ron')
ylabel('radd')
zlabel('roff')
set(gca,'linewidth',1)
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
box on

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure3/Figure3C'])
    print('-dpdf',['./Figures/Figure3/Figure3C'])
end

%% Figure3D

clearvars -except saveFigures
clc

% load('/Volumes/MELANOMA/Data/RareParameters/3nodes/rare_par1000_3_2.mat')
load('./Data/RareParameters/3nodes/rare_par1000_3_2.mat')
R32 = length(rare_par);

% load('/Volumes/MELANOMA/Data/RareParameters/3nodes/Constrained/rare_par1000C_3_2.mat')
load('./Data/RareParameters/3nodes/Constrained/rare_par1000C_3_2.mat')
R32C = length(rare_par);

% load('/Volumes/MELANOMA/Data/RareParameters/5nodes/rare_par1000_5_3.mat')
load('./Data/RareParameters/5nodes/rare_par1000_5_3.mat')
R53 = length(rare_par);

% load('./Data/RareParameters/5nodes/Constrained/rare_par1000C_5_3.mat')
% R53C = length(rare_par);
R53C = 311;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.5;
ylen = 0.5;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
b = bar([R32, R32C; R53, R53C],'FaceColor',[51 51 51]./255, 'EdgeColor','none');
b(2).FaceColor = [255, 127, 42]./255;
xlabel('network')
ylabel('% simulations')
yticks([0, 400])
yticklabels({'0', '40'})
xticklabels({'3.2', '5.3'})
set(gca,'linewidth',2)
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
box off

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure3/Figure3D'])
    print('-dpdf',['./Figures/Figure3/Figure3D'])
end