
addpath(genpath(pwd));

saveFigures = true;

%% Figure2A

clearvars -except saveFigures
clc

%network 3.2, parameter 915

% load('/Volumes/MELANOMA/Data/Simulations/3nodes/S_outpar1000_3_2_10')
% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Simulations/3nodes/S_outpar1000_3_2_10')
load('./Data/Data1000')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;
param = 15;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.25;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(411100-100:412100+100,S_outpar{param}(1:3,411100-100:412100+100),'LineWidth',2);
hold on 
plot(411100-100:412100+100,repmat(thres(915),1,length(411100-100:412100+100)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'},  {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255});
xlim([411100-100, 412100+100])
xticks([])
xticklabels({})
ylim([0, 220])
yticks([0, 220])
yticklabels({'0', '2.2'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
ylabel('gene product (10^2)')
xlabel('time')
box off

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure2/Figure2A'])
    print('-dpdf',['./Figures/Figure2/Figure2A'])
end

%% Figure2B

clearvars -except saveFigures
clc

%left panel 'Data'
%data taken from Schaffer et al., Nature, 2018: Rare cell variability and
%drug-induced reprogramming as a mode of cancer drug resistance, 
%Extended Figure 5 d left panel

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram([zeros(8932,1);ones(720,1);repmat(2,197,1);repmat(3,82,1);...
    repmat(4,37,1);repmat(5,17,1);repmat(6,8,1);repmat(7,5,1);...
    repmat(8,2,1)],'Facecolor',[100,100,100]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1);
set(gca,'linewidth',2)
box off
ylim([0,1])
yticks([0, 0.5, 1])
yticklabels({'0', '50', '100'})
xlim([-0.5,8.5])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
xlabel('no. highly expressed genes')
ylabel('% of cells')
title('data')

%right panel 'Simulation'
% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Data1000')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;

i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
loadsol = sprintf('./Data/CriteriaAnalysis/8nodes/sol10008_1_136');
load(loadsol)
param = 4;
histogram(sol{4}.samp_time,'Facecolor',[0,0,0]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1);
set(gca,'linewidth',2)
box off
ylim([0,1])
yticks([])
xlim([-0.5,8.5])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
title('simulation')

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure2/Figure2B'])
    print('-dpdf',['./Figures/Figure2/Figure2B'])
end

%% Figure2C

clearvars -except saveFigures
clc

%left panel 'Data'
%data taken from Schaffer et al., Nature, 2018: Rare cell variability and
%drug-induced reprogramming as a mode of cancer drug resistance

% data = importfile('/Volumes/MELANOMA/Data/WM9_noDrug_20150618.txt');
data = importfile('./Data/WM9_noDrug_20150618.txt');

%NGFR
pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(data(:,12),'Facecolor',[100,100,100]./255,'Normalization','probability',...
    'Binwidth',50,'EdgeColor','none', 'FaceAlpha', 1);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([0, 0.5, 1])
yticklabels({'0', '50', '100'})
xticks([0, 750])
xlabel('gene product')
ylabel('% of cells')
title('data')
xlim([0,750])
ylim([0 1])

%carpet
i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
plot(data(:,12),ones(length(data(:,12)),1),'.','Color',[100,100,100]./255);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca, 'visible', 'off')
set(0, 'DefaultFigureRenderer', 'painters');
xlim([0,750])

%right panel 'Simulation'
i = 3;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
loadsol = sprintf('./Data/CriteriaAnalysis/8nodes/sol10008_1_136');
load(loadsol)
param = 4;
histogram(sol{4}.samp(1,:),'Facecolor',[135 170 222]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1, 'Binwidth',50);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([])
xticks([0, 750])
xlabel('gene product')
title('simulation')
xlim([0,750])
ylim([0 1])

%carpet
i = 4;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
plot(sol{4}.samp(1,:),ones(1000,1),'.','Color',[135 170 222]./255);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca, 'visible', 'off')
set(0, 'DefaultFigureRenderer', 'painters');
xlim([0,750])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure2/Figure2C'])
    print('-dpdf',['./Figures/Figure2/Figure2C'])
end


%% Figure2D

clearvars -except saveFigures
clc

Net = [2,4,10,80];
ncount = 0;

for n_species = [2,3,5,8]
    ncount = ncount + 1;
    
    R1 = zeros(1,0);
    R2 = zeros(1,0);
    R3 = zeros(1,0);
   
    for inet = 1:Net(ncount)
        rare_par = [];
        
%         loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',...
%             n_species,n_species,inet);
        loadrare = sprintf('./Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',...
            n_species,n_species,inet);
        load(loadrare)
        R1 = [R1,rare_par];
        
%         loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/Replicates2/rare_par10002_%d_%d.mat',...
%             n_species,n_species,inet);
        loadrare = sprintf('./Data/RareParameters/%dnodes/Replicates2/rare_par10002_%d_%d.mat',...
            n_species,n_species,inet);
        load(loadrare)
        R2 = [R2,rare_par];
        
%         loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/Replicates3/rare_par10003_%d_%d.mat',...
%             n_species,n_species,inet);
        loadrare = sprintf('./Data/RareParameters/%dnodes/Replicates3/rare_par10003_%d_%d.mat',...
            n_species,n_species,inet);
        load(loadrare)
        R3 = [R3,rare_par];
    end
    
    RAll(ncount,:) = [length(R1),length(R2),length(R3)];
end

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.25;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
boxplot([RAll(1,:)',RAll(2,:)', RAll(3,:)', RAll(4,:)'],'Colors', [0,0,0]./255)
ylim([0,240])
xlabel('network size')
xticklabels({'2','3','5','8'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
ylabel('simulations')
yticks([0,120,240])
set(gca,'linewidth',2)
box off
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure2/Figure2D'])
    print('-dpdf',['./Figures/Figure2/Figure2D'])
end

%% Figure2E

clearvars -except saveFigures
clc

n_species = 5;
for inet = 1:10
    rare_par = [];
%     loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',...
%                 n_species,n_species,inet);
    loadrare = sprintf('./Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',...
        n_species,n_species,inet);
    load(loadrare)
    R(inet) = length(rare_par);
end

% load('/Volumes/MELANOMA/Data/M_iso5')
load('./Data/M_iso5')
for jnet = 1:length(M_iso)
    C(jnet) = sum(M_iso{jnet}(:,1));
end

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.25;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
plot(C,R,'.','Markersize',20,'Color', [211,95,95]./255)
ylabel('simulations')
xlabel('connectivity')
xticks([1,2,3,4,5])
yticks([0,20,40,60])
xlim([0.5,5.5])
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca,'linewidth',2)

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure2/Figure2E'])
    print('-dpdf',['./Figures/Figure2/Figure2E'])
end

%% Figure2F

count = 1;
for n_species = [2,3,5,8]
    
    loadmat = sprintf('./Data/M_iso%d',n_species);
%     loadmat = sprintf('/Volumes/MELANOMA/Data/M_iso%d',n_species);
    load(loadmat);
    
    for i = 1:length(M_iso)
        
        conn{count}(i) = sum(M_iso{i}(1,:));
        
        loops{count}(i) =  M_iso{i}(1,1);
        
%         loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',n_species,n_species,i);
        loadrare = sprintf('./Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',n_species,n_species,i);
        load(loadrare)
        
        Rare_par{i} = rare_par;
        
        rare{count}(i) = length(rare_par);      
    end
    
    R_noloops{count} = rare{count}(loops{count} == 0);
    R_loops{count} = rare{count}(loops{count} == 1);
    
     Ratio{count} = rare{count}(loops{count} == 0)./rare{count}(loops{count} == 1);
    
    count = count + 1;
end
A = cell2mat(R_loops);
B = cell2mat(R_noloops);
C = cell2mat(Ratio);
%plot of network size

jitter = 1/100.*randn(length(C(C>0 & C < 100)),1) + 1;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.4;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
[h,L,MX,MED]=violin(C(C>0 & C < 100)', 'facecolor', [50,50,50;50,50,50]./255);
hold on
plot(jitter,C(C>0 & C < 100),'.','Color','k','Markersize',15)
hold on
plot([0,2],[1,1],':','Color','k','LineWidth',2)
xticks([])
yticks([0:6])
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
box off
set(gca,'linewidth',2)
ylabel('fold change')

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/Figure2/Figure2F'])
    print('-dpdf',['./Figures/Figure2/Figure2F'])
end
%% Figure2G

clearvars -except saveFigures
clc

%0 - stabel low
%1 - uncorrelated transient high
%2 - rare transient correlaetd high
%3 - stably high

Net = [2,4,10];
netcount = 0;
C = [];

for n_species = [2,3,5]
    netcount = netcount +1;
    loadclasses = sprintf('./Data/ClassAnalysis/classAll%d',n_species);
    load(loadclasses)
    
    ind = find(ismember(class(:,3),[48,52,95,152,282,408,558,769,915,994]) == 1);
    
    C1 = class(ind,4);
    
    C(1:10,size(C,2)+1:sum(Net(1:netcount))) = reshape(C1,[10,(Net(netcount))]);
end

flipud(C)

