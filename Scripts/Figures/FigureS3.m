
addpath(genpath(pwd));

saveFigures = true;

%% SupplementaryFigure3A

clearvars -except saveFigures
clc

% data = importfile('/Volumes/MELANOMA/Data/WM9_noDrug_20150618.txt');
data = importfile('./Data/WM9_noDrug_20150618.txt');

lim = [35, 100, 400, 750, 100, 140, 300, 85, 190];

count = 1;

for i = [4,8,11,12,13,15,21,18,16]
   
    pos = [0.1,0.32,0.54,0.76];
    xlen = 0.2;
    ylen = 0.2;
    
    figure
    j = 1;
    subplot('Position',[pos(j),pos(1),xlen,ylen]);
    histogram(data(:,i),'Facecolor',[100,100,100]./255,'Normalization','probability',...
        'Binwidth', lim(count)/30,'EdgeColor','none', 'FaceAlpha', 1);
    set(gca,'linewidth',2)
    box off
    set(0,'DefaultAxesFontName','Arial');
    set(gca,'FontSize',12)
    yticks([0, 0.5, 1])
    yticklabels({'0', '50', '100'})
    xlabel('gene product')
    ylabel('% of cells')
    ylim([0 1])
    xlim([0,lim(count)])
    
    %carpet
    j = 2;
    subplot('Position',[pos(j),pos(1),xlen,ylen]);
    plot(data(:,i),ones(length(data(:,12)),1),'.','Color',[100,100,100]./255);
    set(gca,'linewidth',2)
    box off
    set(0,'DefaultAxesFontName','Arial');
    set(gca,'FontSize',12)
    set(gca, 'visible', 'off')
    set(0, 'DefaultFigureRenderer', 'painters');
    xlim([0,lim(count)])
    
    if saveFigures
        set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%         path = sprintf('/Volumes/MELANOMAII/Figures/FigureS3/FigureS3A%d',i);
        path = sprintf('./Figures/FigureS3/FigureS3A%d',i);
        print('-dpdf',[path])
    end

    count = count + 1;
end

%% SupplementaryFigure3B

clearvars -except saveFigures
clc

n_species = 3;
Net = 1;
Subnet = 10;

count = 1;
countfig = 0;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure


for inet = [1,1]
    
    countfig = countfig + 1;
%     loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/Replicates3/rare_par10003_%d_%d.mat',...
%         n_species,n_species,inet);
    loadrare = sprintf('./Data/RareParameters/%dnodes/Replicates3/rare_par10003_%d_%d.mat',...
        n_species,n_species,inet);
    load(loadrare)
    
    Subnet = [1,3];
    
    for isubnet = Subnet(countfig)
        
%         r1 = rare_par(rare_par > isubnet*1000/Subnet-1000/Subnet);
%         r2 = rare_par(rare_par <= isubnet*1000/Subnet);
%         rare_par_isubnet = intersect(r1,r2);
        
%         if isempty(rare_par_isubnet) == 0
            
            loadsol = sprintf('./Data/CriteriaAnalysis/%dnodes/Replicates3/sol10003%d_%d_%d.mat',n_species,n_species,inet,isubnet);
            load(loadsol)
            
            R = [27,239];
            
            for i = R(countfig)
                
                DataSim = sol{i-(length(sol)*isubnet-length(sol))}.samp(1,:);
                PrctlSim = prctile(sol{i-(length(sol)*isubnet-length(sol))}.samp(1,:),99);
    
                a = histcounts(DataSim,'BinWidth',1);
                bin = 1;
                
                if find(a == max(a)) <= 15

                    anew = max(a);
                    error = 100;
                    countwhile = 0;
                    while abs(error) > 1
                        if error > 0
                            anew = anew + 10;
                        else
                            anew = anew - 10;
                            if anew < 0
                                DataExp = zeros(0,1);
                                break
                            end
                        end
                        DataExp = exprnd(anew,1,1000);
                        b = histcounts(DataExp,'BinWidth',bin);
                        error = max(b)-max(a);
                        countwhile = countwhile + 1;
                        if countwhile > 1000
                            DataExp = zeros(0,1);
                            break
                        end
                    end

                    A(count) = anew;
                    if isempty(DataExp) == 0
                        PrctlExp = prctile(DataExp,99);
                        
                        j = 1;
                        subplot('Position',[pos(j),pos(1),xlen,ylen]);
                        
                        histogram(DataSim,'Binwidth',10,'Normalization','probability',...
                            'EdgeColor','none','Facecolor', [211,95,95]./255, 'FaceAlpha', 200/255)
                        hold on
                        histogram(DataExp,'Binwidth',10,'Normalization','probability',...
                            'EdgeColor','none','Facecolor', [150,150,150]./255, 'FaceAlpha', 150/255)
                        box off
                        set(gca,'linewidth',2)
                        xlabel('gene product')
                        ylabel('% of cells')
                        xticks([0,700])
                        yticks([0,0.14])
                        yticklabels({'0','14'})
                        xlim([0,700])
                        ylim([0,0.14])
                        set(0,'DefaultAxesFontName','Arial'); 
                        set(gca,'FontSize',12)

                        Com(count,:) = [1,n_species,inet,isubnet,i-(length(sol)*isubnet-length(sol)),i,PrctlSim,PrctlExp];

                        count = count+1;

                    else
                        Com(count,:) = [1,n_species,inet,isubnet,i-(length(sol)*isubnet-length(sol)),i,NaN,NaN]; %0 for not run through
                        count = count+1;
                    end
                else
                    Com(count,:) = [0,n_species,inet,isubnet,i-(length(sol)*isubnet-length(sol)),i,NaN,NaN];  %NaN for max bin > 15
                    count = count+1;

                    anew = max(a);
                    error = 100;
                    countwhile = 0;
                    while abs(error) > 1
                        if error > 0
                            anew = anew + 10;
                        else
                            anew = anew - 10;
                            if anew < 0
                                DataExp = zeros(0,1);
                                break
                            end
                        end
                        DataExp = exprnd(anew,1,1000);

                        b = histcounts(DataExp,'BinWidth',bin);
                        error = max(b)-max(a);
                        countwhile = countwhile + 1;
                        if countwhile > 1000
                            DataExp = zeros(0,1);
                            break
                        end
                    end
                    
                    j = 2;
                    subplot('Position',[pos(j),pos(1),xlen,ylen]);
                    
                    histogram(sol{i-(length(sol)*isubnet-length(sol))}.samp(1,:),...
                        'Binwidth',50,'Normalization','probability','EdgeColor','none',...
                        'Facecolor', [211,95,95]./255, 'FaceAlpha', 200/255)
                    hold on
                    histogram(DataExp,'Binwidth',50,'Normalization','probability',...
                        'EdgeColor','none','Facecolor', [150,150,150]./255,...
                        'FaceAlpha', 150/255)
                    box off
                    set(gca,'linewidth',2)
                    xlabel('gene product')
                    xticks([0,2200])
                    yticks([])
                    xlim([0,2200])
                    ylim([0,0.14])
                    set(0,'DefaultAxesFontName','Arial');
                    set(gca,'FontSize',12)
                    
                end

            end
%         end
    end
end

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS3/FigureS3B'])
    print('-dpdf',['./Figures/FigureS3/FigureS3B'])
end

%% SupplementaryFigure3C

clearvars -except saveFigures
clc

%left panel

% data = importfile('/Volumes/MELANOMA/Data/WM9_noDrug_20150618.txt');
data = importfile('./Data/WM9_noDrug_20150618.txt');

count = 1;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.5;
ylen = 0.5;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);

for i = [4,8,11,12,13,15,21,18,16]
    
    DataSim = data(:,i);
    PrctlSim = prctile(data(:,i),99);
    
    a = histcounts(DataSim,'BinWidth',1);
    bin = 1;
    
    anew = mean(DataSim);
    error = 100;
    countwhile = 0;
    while abs(error) > 1
        if error > 0
            anew = anew + 0.1;
        else
            anew = anew - 0.1;
            if anew < 0
                DataExp = zeros(0,1);
                break
            end
        end
        DataExp = exprnd(anew,1,length(DataSim));    
        b = histcounts(DataExp,'BinWidth',bin);
        error = max(b)-max(a);
        countwhile = countwhile + 1;
        if countwhile > 100000
            DataExp = zeros(0,1);
            break
        end
    end

    A(count) = anew;
    
%     figure
%     histogram(DataSim,'Binwidth',1)
%     hold on
%     histogram(DataExp,'Binwidth',1)
    
    if isempty(DataExp) == 0 
        PrctlExp = prctile(DataExp,99);
        
        Com(count,:) = [i,PrctlSim,PrctlExp];
        count = count+1;
        
    else
        Com(count,:) = [i,NaN,NaN]; %0 for not run through
        count = count+1;
    end
end

plot(Com(:,3),Com(:,2),'.','MarkerSize',10, 'Color', [211,95,95]./255)
hold on
plot([1,300],[1,300],':','Color','k','LineWidth',2)
% box off
% set(gca,'linewidth',2)
% xlabel('99th percentile of fitted exponential distribution')
% ylabel({'99th percentile of expression'; 'distributions of simulations with'; 'rare coordinated high states'})
% xticks([0,300])
% yticks([0,300])
% set(0,'DefaultAxesFontName','Arial');
% set(gca,'FontSize',12)
% xlim([0,300])
% ylim([0,300])

load('./Data/Com1.mat')
plot(Com(:,8),Com(:,7),'.','MarkerSize',10, 'Color', [200,200,200]./255)
hold on
plot([1,300],[1,300],':','Color','k','LineWidth',2)
box off
set(gca,'linewidth',2)
xlabel('99th percentile of fitted exponential distribution')
ylabel({'99th percentile of expression distributions'; 'of experimental data and simulations'})
xticks([0,300])
yticks([0,300])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
xlim([0,300])
ylim([0,300])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS3/FigureS3Cleft'])
    print('-dpdf',['./Figures/FigureS3/FigureS3Cleft'])
end

%right panel
pos = [0.1,0.32,0.54,0.76];
xlen = 0.5;
ylen = 0.5;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
plot(Com(:,8),Com(:,7),'.','MarkerSize',10, 'Color', [200,200,200]./255)
hold on
plot([1,3000],[1,3000],':','Color','k','LineWidth',2)
box off
set(gca,'linewidth',2)
xlabel('99th percentile of fitted exponential distribution')
xticks([0,3000])
yticks([0,3000])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
xlim([0,3000])
ylim([0,3000])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS3/FigureS3Cright'])
    print('-dpdf',['./Figures/FigureS3/FigureS3Cright'])
end

%% SupplementaryFigure3D

%analysis

% Net = [2,4,10,80];
% Subnet = [10,10,100,250];
% Nspecies = [2,3,5,8];
% 
% countgini = 1;
% countginiII = 1;
% count = 1;
% 
% load('/Volumes/MELANOMA/Data/Data1000')
% thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;
% 
% for n_species = [2,3,5,8]
%     % for n_species = 8
%     n_species
%     for inet = 1:Net(count)
%         inet
%         
%         loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',n_species,n_species,inet);
%         load(loadrare)
%         
%         for isubnet = 1:Subnet(count)
%             isubnet
%             r1 = rare_par(rare_par > isubnet*1000/Subnet(count)-1000/Subnet(count));
%             r2 = rare_par(rare_par <= isubnet*1000/Subnet(count));
%             rare_par_isubnet = intersect(r1,r2);
%             
%             if isempty(rare_par_isubnet) == 0
%                 
%                 loadsol = sprintf('/Volumes/MELANOMAII/Data/CriteriaAnalysis/%dnodes/sol1000%d_%d_%d.mat',n_species,n_species,inet,isubnet);
%                 load(loadsol)
%                 
%                 for param = rare_par_isubnet
%                     
%                     genecount = sol{param-(isubnet-1)*1000/Subnet(count)}.samp(1,:);
%                     
%                     G(countgini) = gini(ones(1000,1),genecount);
%                     countgini = countgini + 1;
%                 end
%             end
%             
%             if isempty(rare_par_isubnet) == 1
%                 
%                 loadsol = sprintf('/Volumes/MELANOMAII/Data/CriteriaAnalysis/%dnodes/sol1000%d_%d_%d.mat',n_species,n_species,inet,isubnet);
%                 load(loadsol)
%                 
%             end
%             
%             for paramII = setdiff(1000/Subnet(count)*(isubnet-1)+1:1000/Subnet(count)*isubnet,rare_par_isubnet)
%                 
%                 genecount = sol{paramII-(isubnet-1)*1000/Subnet(count)}.samp(1,:);
%                 
%                 GII(countginiII) = gini(ones(1000,1),genecount);
%                 countginiII = countginiII + 1;
%             end
%             
%             
%         end
%     end
%     count = count + 1;
% end

load('./Data/G')
load('./Data/GII')

% data = importfile('/Volumes/MELANOMA/Data/WM9_noDrug_20150618.txt');
data = importfile('./Data/WM9_noDrug_20150618.txt');

count = 1;
for i = [4,8,11,12,13,15,21,18,16]
    
    GExp(count) = gini(ones(length(data(:,i)),1),data(:,i));
    count = count + 1;
end

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
histogram(G,'Binwidth',0.025,'Normalization','probability','Facecolor',[211,95,95]./255,'Edgecolor','none','FaceAlpha',1)
hold on
histogram(GII,'Binwidth',0.025,'Normalization','probability','Facecolor',[100,100,100]./255,'Edgecolor','none','FaceAlpha',1)
hold on
for j = 1:9
    plot([GExp(j),GExp(j)],[0,0.05],'Color','k','Linewidth',1)
    hold on
end
xticks([0,0.5,1])
yticks([0,0.2,0.4])
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
box off
set(gca,'linewidth',2)
ylabel('probability')
xlabel('Gini coefficient')
ylim([0,0.4])
xlim([0,1])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS3/FigureS3D'])
    print('-dpdf',['./Figures/FigureS3/FigureS3D'])
end


%% SupplementaryFigure3E

%network 2.1
% load('/Volumes/MELANOMAII/Data/QuantificationAnalysis/2nodes/solQuant10002_1')
load('./Data/QuantificationAnalysis/2nodes/solQuant10002_1')
s21 = solQuant.AllTimeJack([5,11,16,22,34,59,60,64])/max(solQuant.AllTimeJack([5,11,16,22,34,59,60,64]));

%network 2.2
% load('/Volumes/MELANOMAII/Data/QuantificationAnalysis/2nodes/solQuant10002_2')
load('./Data/QuantificationAnalysis/2nodes/solQuant10002_2')
s22 = solQuant.AllTimeJack([4,8,12,15,20,33,34,36])/max(solQuant.AllTimeJack([4,8,12,15,20,33,34,36]));

%network 3.1
% load('/Volumes/MELANOMAII/Data/QuantificationAnalysis/3nodes/solQuant10003_1')
load('./Data/QuantificationAnalysis/3nodes/solQuant10003_1')
s31 = solQuant.AllTimeJack/max(solQuant.AllTimeJack);

%network 3.2
% load('/Volumes/MELANOMAII/Data/QuantificationAnalysis/3nodes/solQuant10003_2')
load('./Data/QuantificationAnalysis/3nodes/solQuant10003_2')
s32 = solQuant.AllTimeJack/max(solQuant.AllTimeJack);

%network 3.3
% load('/Volumes/MELANOMAII/Data/QuantificationAnalysis/3nodes/solQuant10003_3')
load('./Data/QuantificationAnalysis/3nodes/solQuant10003_3')
s33 = solQuant.AllTimeJack/max(solQuant.AllTimeJack);

%network 5.2
% load('/Volumes/MELANOMAII/Data/QuantificationAnalysis/5nodes/solQuant10005_2')
load('./Data/QuantificationAnalysis/5nodes/solQuant10005_2')
s52 = solQuant.AllTimeJack/max(solQuant.AllTimeJack);

%network 5.3
% load('/Volumes/MELANOMAII/Data/QuantificationAnalysis/5nodes/solQuant10005_3')
load('./Data/QuantificationAnalysis/5nodes/solQuant10005_3')
s53 = solQuant.AllTimeJack/max(solQuant.AllTimeJack);

M = [s21',s22',s31',s32',s33',s52',s53'];

for i = 1:size(M,2)
    
    Max = max(M(:,i));
    for j = 1:size(M,1)
        count = M(j,i);
        red(j,i) = 160 +  abs(count-Max)*(244-160);
        green(j,i) = 44 + abs(count-Max)*(215-44);
        blue(j,i) = 44 + abs(count-Max)*(215-44);
    end
end

%RGB values
red
green
blue

