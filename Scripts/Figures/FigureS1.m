
addpath(genpath(pwd));

saveFigures = true;

%% SupplementaryFigure1A

clearvars -except saveFigures
clc

load('./Data/Simulations/3nodes/S_outparModelOld3_1_1')
load('./Data/Simulations/3nodes/R_outparModelOld3_1')

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);

%SupplementaryFigure1C left

param = 11;

h = plot(54000:68000,S_outpar{param}(1:3,54000:68000),'LineWidth',2);
set(h, {'color'}, {[197 90 17]./255; [68 114 196]./255; [255, 192 0]./255});
hold on 
plot(54000:68000,repmat(R(param,10)/0.9*0.8,1,length(54000:68000)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'}, {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255});
set(gca,'linewidth',2)
box off
xticks([54000, 68000])
yticks([0, 40, 80])
xticklabels({'0', '14000'})
yticks([0, 40, 80])
ylim([0,80])
xlim([54000, 68000])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
ylabel('gene product')
xlabel('time')
box off

%SupplementaryFigure1C middle

param = 48;

%plot simulation
i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(54000:68000,S_outpar{param}(1:3,54000:68000),'LineWidth',2);
set(h, {'color'}, {[197 90 17]./255; [68 114 196]./255; [255, 192 0]./255});
hold on 
plot(54000:68000,repmat(R(param,10)/0.9*0.8,1,length(54000:68000)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'}, {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255});
set(gca,'linewidth',2)
xlim([54000, 68000])
box off
xticks([54000, 68000])
yticks([])
xticklabels({'0', '14000'})
yticks([])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
xlabel('time')
box off

%SupplementaryFigure1C right

param = 3;

i = 3;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
h = plot(54000:68000,S_outpar{param}(1:3,54000:68000),'LineWidth',2);
set(h, {'color'}, {[197 90 17]./255; [68 114 196]./255; [255, 192 0]./255});
hold on 
plot(54000:68000,repmat(R(param,10)/0.9*0.8,1,length(54000:68000)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'}, {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255});
set(gca,'linewidth',2)
box off
xlim([54000, 68000])
xticks([54000, 68000])
yticks([])
xticklabels({'0', '14000'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
xlabel('time')
box off

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS1/FigureS1A'])
    print('-dpdf',['./Figures/FigureS1/FigureS1A'])
end

%% SupplementaryFigure1B

clearvars -except saveFigures
clc

load('./Data/Simulations/3nodes/S_outparModelOld3_1_1')
load('./Data/Simulations/3nodes/R_outparModelOld3_1')
% load('./Data/Simulations/3nodes/S_outparModelOld3_1_1')
% load('./Data/Simulations/3nodes/R_outparModelOld3_1')

n_species = 3;
inet = 1;
doplot = 'no';
tcell = 1000;                        
param = 48;

isubnet = 1;
for i = param
    param = (isubnet-1)*length(S_outpar)+i;
    threshold = R(i,10)/0.9*0.8;
    [maxjackpot_sol,desc_sol,rightskewed_sol,unimodal_sol,samp_sol,samp_time_sol,rand_time_sol] ...
        = analyzeQual(n_species,tcell,threshold,S_outpar{i});
    sol{i}.maxjackpot = maxjackpot_sol;
    sol{i}.desc = desc_sol;
    sol{i}.rightskewed = rightskewed_sol;
    sol{i}.unimodal = unimodal_sol;
    sol{i}.samp = samp_sol;
    sol{i}.samp_time = samp_time_sol;
    sol{i}.time = rand_time_sol;
end

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{param}.samp_time,'Normalization','probability','Edgecolor','none',...
    'Facecolor',[0,0,0],'Facealpha',1)
xlim([-0.5,3.5])
ylim([0,1])
yticks([0,1])
yticklabels({'0','100'})
xlabel('no. highly expressed genes')
xticks([0,1,2,3])
ylabel('% of cells')
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)

i = 2;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
histogram(sol{48}.samp,'Normalization','probability','Edgecolor','none',...
    'Facecolor',[135 170 222]./255,'Facealpha',1,'Binwidth',3)
xlabel('gene product')
xlim([0,90])
xticks([0,90])
ylim([0,1])
yticks([])
yticklabels({'0','100'})
alpha(1)
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS1/FigureS1B'])
    print('-dpdf',['./Figures/FigureS1/FigureS1B'])
end

%% SupplementaryFigure1C

clearvars -except saveFigures
clc

for i = 2:10
    %    nameload = sprintf('/Volumes/MELANOMA/Data/M_iso%d',i);
    nameload = sprintf('./Data/M_iso%d',i);
    load(nameload);
    L(i) = log10(length(M_iso));
    
    N(i) = log10(2^(i*i));
end
N(1) = log10(2);
L(1) = log10(2);

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
plot(1:10,N,'-','Linewidth',2,'Color',[150,150,150]./255)
hold on
plot(1:10,N,'.','Markersize',12,'Color',[150,150,150]./255)
hold on
plot(1:10,L,'-','LineWidth',2,'Color', 'k')
hold on
plot(1:10,L,'.','Markersize',12,'Color','k')
set(gca,'linewidth',2)
box off
xlabel('network size')
ylabel('log_{10}(networks)')
xlim([1,10])
ylim([0,40])
yticks([0,20,40])
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS1/FigureS1C'])
    print('-dpdf',['./Figures/FigureS1/FigureS1C'])
end

%% SupplementaryFigure1D

clearvars -except saveFigures
clc

importFreq1;
Freq1 = WM989noDrug20150618Summary;
importFreq2;
Freq2 = WM989noDrug20150618Summary;

for i = 1:size(Freq2,1)
    for j = 1:size(Freq2,2)
        
        if isnan(Freq2(i,j)) == 0
%             Freq2norm(i,j) = Freq2(i,j)/(Freq1(i)*Freq1(j));
              Freq2norm(i,j) = Freq2(i,j);
        else
            Freq2norm(i,j) = NaN;
        end
        
    end
end

%remove NaN columns
Freq2normNaN = [Freq2norm([1:3,5:7,9:11,13],[1:3,5:7,9:11,13,14]);nan(1,11)];

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.4;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
b = imagesc(Freq2normNaN);
set(b,'AlphaData',~isnan(Freq2normNaN))
colormap(flipud(gray))
set(gca,'linewidth',2)
yticks([])
xticks([])
set(gca,'linewidth',2)
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
colorbar

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS1/FigureS1D'])
    print('-dpdf',['./Figures/FigureS1/FigureS1D'])
end

%% SupplementaryFigure1E

%Hierarchical

% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Data1000')

thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;

Subnet = 10;
count = 1;

for isubnet = 1:Subnet
   
    loadrare = sprintf('./Data/Hierarchical/rare_par1000Hierarchical_5_1_%d.mat',isubnet);
    load(loadrare)
    
    r1 = rare_par(rare_par > isubnet*1000/Subnet(count)-1000/Subnet(count));
    r2 = rare_par(rare_par <= isubnet*1000/Subnet(count));
    rare_par_isubnet = intersect(r1,r2);
    
    if isempty(rare_par_isubnet) == 0
        
        loadsol = sprintf('./Data/Hierarchical/sol1000Hierarchical_5_1_%d',isubnet);
        load(loadsol)
        
        for irare = rare_par_isubnet
            
            irare
            
            for icell = 1:5

                count1 = 0;
                count2 = 0;
                count3 = 0;
                count4 = 0;
                count5 = 0;
                
                FreqSingle{irare}(icell) = sum(sol{irare-(isubnet-1)*100}.samp(icell,:) > thres(irare))/1000;
                
                PosSingle = find(sol{irare-(isubnet-1)*100}.samp(icell,:) > thres(irare));
                
                for ipos = PosSingle
                    
                    Double = find(sol{irare-(isubnet-1)*100}.samp(:,ipos)  > thres(irare));
                    if length(Double) > 1   
                        if ismember(1,Double) == 1
                            count1 = count1 + 1;
                        end
                        if ismember(2,Double) == 1
                            count2 = count2 + 1;
                        end
                        if ismember(3,Double) == 1
                            count3 = count3 + 1;
                        end
                        if ismember(4,Double) == 1
                            count4 = count4 + 1;
                        end
                        if ismember(5,Double) == 1
                            count5 = count5 + 1;
                        end
                    end
                    
                end
                
                FreqDouble{irare}(icell).cell1 = count1/1000;
                FreqDouble{irare}(icell).cell2 = count2/1000;
                FreqDouble{irare}(icell).cell3 = count3/1000;
                FreqDouble{irare}(icell).cell4 = count4/1000;
                FreqDouble{irare}(icell).cell5 = count5/1000;
            end
            
            
        end
        
    end
       
end

%remove NaN columns
Freq2 = [NaN, FreqDouble{14}(1).cell2, FreqDouble{14}(1).cell3, FreqDouble{14}(1).cell4, FreqDouble{14}(1).cell5;...
    NaN, NaN, FreqDouble{14}(2).cell3, FreqDouble{14}(2).cell4, FreqDouble{14}(2).cell5;...
    NaN, NaN, NaN, FreqDouble{14}(3).cell4, FreqDouble{14}(3).cell5;...
    NaN, NaN, NaN, NaN FreqDouble{14}(4).cell5;...
    NaN, NaN, NaN, NaN, NaN];

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.4;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
b = imagesc(Freq2);
set(b,'AlphaData',~isnan(Freq2))
colormap(flipud(gray))
set(gca,'linewidth',2)
yticks([])
xticks([])
set(gca,'linewidth',2)
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
colorbar

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS1/FigureS1E'])
    print('-dpdf',['./Figures/FigureS1/FigureS1E'])
end

%% SupplementaryFigure1F

isubnet = 97;

% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Data1000')

thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;


% load('/Volumes/MELANOMA/Data/RareParameters/5nodes/rare_par1000_5_3.mat');
load('./Data/RareParameters/5nodes/rare_par1000_5_3.mat');

load('./Data/CriteriaAnalysis/5nodes/sol10005_1_97');


for irare = 968
    
    for icell = 1:5
        
        count1 = 0;
        count2 = 0;
        count3 = 0;
        count4 = 0;
        count5 = 0;
        
        FreqSingle{irare}(icell) = sum(sol{irare-(isubnet-1)*10}.samp(icell,:) > thres(irare))/1000;
        
        PosSingle = find(sol{irare-(isubnet-1)*10}.samp(icell,:) > thres(irare));
        
        for ipos = PosSingle
            
            Double = find(sol{irare-(isubnet-1)*10}.samp(:,ipos)  > thres(irare));
            if length(Double) > 1
                if ismember(1,Double) == 1
                    count1 = count1 + 1;
                end
                if ismember(2,Double) == 1
                    count2 = count2 + 1;
                end
                if ismember(3,Double) == 1
                    count3 = count3 + 1;
                end
                if ismember(4,Double) == 1
                    count4 = count4 + 1;
                end
                if ismember(5,Double) == 1
                    count5 = count5 + 1;
                end
            end
            
        end
        
        FreqDouble{irare}(icell).cell1 = count1/1000;
        FreqDouble{irare}(icell).cell2 = count2/1000;
        FreqDouble{irare}(icell).cell3 = count3/1000;
        FreqDouble{irare}(icell).cell4 = count4/1000;
        FreqDouble{irare}(icell).cell5 = count5/1000;
    end
    
    
end
        
       

%remove NaN columns
Freq2 = [NaN, FreqDouble{968}(1).cell2, FreqDouble{968}(1).cell3, FreqDouble{968}(1).cell4, FreqDouble{968}(1).cell5;...
    NaN, NaN, FreqDouble{968}(2).cell3, FreqDouble{968}(2).cell4, FreqDouble{968}(2).cell5;...
    NaN, NaN, NaN, FreqDouble{968}(3).cell4, FreqDouble{968}(3).cell5;...
    NaN, NaN, NaN, NaN FreqDouble{968}(4).cell5;...
    NaN, NaN, NaN, NaN, NaN];

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.4;

figure
i = 1;
subplot('Position',[pos(i),pos(1),xlen,ylen]);
b = imagesc(Freq2);
set(b,'AlphaData',~isnan(Freq2))
colormap(flipud(gray))
set(gca,'linewidth',2)
yticks([])
xticks([])
set(gca,'linewidth',2)
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
colorbar

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS1/FigureS1F'])
    print('-dpdf',['./Figures/FigureS1/FigureS1F'])
end
