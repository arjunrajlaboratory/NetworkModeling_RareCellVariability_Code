
addpath(genpath(pwd));

saveFigures = true;

%% SupplementeryFigure6A

clearvars -except saveFigures
clc

Net = [2,4,10,80];

count = 0;
pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure

for n_species = [2,3,5,8]
    
    clearvars -except Net n_species count saveFigures pos xlen ylen
    count = count + 1;
    
    net = Net(count);

    R = zeros(1,0);
    for inet = 1:net
    
%         loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',...
%             n_species,n_species,inet);
        loadrare = sprintf('./Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',...
            n_species,n_species,inet);
        load(loadrare)
        
        R = [R,rare_par];
    
    end

    [hist1,bin] = histcounts(R,'Binwidth',1);

    Z = zeros(1,1000);
    Z(unique(sort(R))) = hist1(hist1>0);
    
    j = count;
    subplot('Position',[pos(j),pos(1),xlen,ylen]);
    plot(1:1000,Z/net,'.','Color',[50,50,50]./255,'Markersize',10)
    hold on
    plot([92,915,183,26,133,544,702,968],Z([92,915,183,26,133,544,702,968])/net,'.',...
        'Color',[255,127,42]./255,'Markersize',10)
    ylim([0,1])
    xticks([0,1000])
    if j == 1
        yticks([0,1])
        yticklabels({'0','100'})
        ylabel({'% of networks'})
    else
        yticks([])
    end
    if j == 2
        xlabel('parameter sets')
    else
        yticks([])
    end
    set(0,'DefaultAxesFontName','Arial'); 
    set(gca,'FontSize',12,'linewidth',2)
    box off
end

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS6/FigureS6A'])
    print('-dpdf',['./Figures/FigureS6/FigureS6A'])
end

%% SupplementeryFigure6C

%for phase space: function PhaseSpace

% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Data1000')
I = [26 92 133 183 544 702 915 968];

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
plot(Data1000(:,3),Data1000(:,5),'.','MarkerSize',5,'Color',[51, 51, 51]./255)
hold on
plot(Data1000(I,3),Data1000(I,5),'.','MarkerSize',10,'Color',[255,127,42]./255)
xlim([0,0.1])
ylim([0,1])
xticks([0,0.1])
yticks([0,1])
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
xlabel({'r_{on}'})
ylabel('r_{add}')
set(gca,'linewidth',2)
box off

j = 2;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
plot(Data1000(:,3),Data1000(:,6),'.','MarkerSize',5,'Color',[51, 51, 51]./255)
hold on
plot(Data1000(I,3),Data1000(I,6),'.','MarkerSize',10,'Color',[255,127,42]./255)
xlim([0,0.1])
ylim([0,0.1])
xticks([0,0.1])
yticks([0,0.1])
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
xlabel({'r_{on}'})
ylabel('r_{off}')
set(gca,'linewidth',2)
box off

j = 3;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
plot(Data1000(:,6),Data1000(:,5),'.','MarkerSize',5,'Color',[51, 51, 51]./255)
hold on
plot(Data1000(I,6),Data1000(I,5),'.','MarkerSize',10,'Color',[255,127,42]./255)
xlim([0,0.1])
ylim([0,1])
xticks([0,0.1])
yticks([0,1])
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
xlabel({'r_{off}'})
ylabel('r_{add}')
set(gca,'linewidth',2)
box off

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS6/FigureS6C'])
    print('-dpdf',['./Figures/FigureS6/FigureS6C'])
end

%% SupplementeryFigure6D

% load('/Volumes/MELANOMA/Data/Data1000')
load('./Data/Data1000')
I = [26 92 133 183 544 702 915 968];
DataSens = repmat(Data1000(I(1),:),70,1);
DataSens(1:10,1) = linspace(0.01,1,10);
DataSens(11:20,2) = linspace(0.001,0.1,10);
DataSens(21:30,3) = linspace(0.001,0.1,10);
DataSens(31:40,4) = linspace(0.1,10,10);
DataSens(41:50,5) = linspace(0.1,1,10);
DataSens(51:60,6) = linspace(0.01,0.1,10);
DataSens(61:70,8) = linspace(2,100,10);
DataSens(:,7) = 0.95*DataSens(:,1)./DataSens(:,2).*DataSens(:,8);

R = [];
for ipar = 1:length(I)
    rareload = sprintf('rare_par1000SensAna%d_3_3_2',I(ipar));
    load(rareload);
    R = [R, rare_par];
end

for jpar = 1:70
    L(jpar) = length(find(R == jpar));
end
L = L/length(I);

for iplot = 1:7

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
if iplot < 7
    plot(DataSens((iplot-1)*10+1:(iplot-1)*10+10,iplot),L((iplot-1)*10+1:(iplot-1)*10+10),'.','Markersize',10,'Color','k')
    hold on
    plot(DataSens((iplot-1)*10+1:(iplot-1)*10+10,iplot),L((iplot-1)*10+1:(iplot-1)*10+10),'-','Linewidth',1,'Color','k')
    hold on
    for j_rarepar = I
        plot([Data1000(j_rarepar,iplot),Data1000(j_rarepar,iplot)],[0,1],'-','Linewidth',1,'Color',[255 130 46]./255)
        hold on
    end
    xticks([0,max(DataSens((iplot-1)*10+1:(iplot-1)*10+10,iplot))])
else
    plot(DataSens((iplot-1)*10+1:(iplot-1)*10+10,8),L((iplot-1)*10+1:(iplot-1)*10+10),'.','Markersize',10,'Color','k')
    hold on
    plot(DataSens((iplot-1)*10+1:(iplot-1)*10+10,8),L((iplot-1)*10+1:(iplot-1)*10+10),'-','Linewidth',1,'Color','k')
    hold on
    for j_rarepar = I
        plot([Data1000(j_rarepar,8),Data1000(j_rarepar,8)],[0,1],'-','Linewidth',1,'Color',[255 130 46]./255)
        hold on
    end
    xticks([0,max(DataSens((iplot-1)*10+1:(iplot-1)*10+10,8))])
end
yticks([0,1])
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
box off
set(gca,'linewidth',2)
ylabel('frequency')
xlabel('parameter')
ylim([0,1])


if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     path = sprintf('/Volumes/MELANOMAII/Figures/FigureS6/FigureS6D%d',iplot);
    path = sprintf('./Figures/FigureS6/FigureS6D%d',iplot);
    print('-dpdf',[path])
end

end

%% SupplementeryFigure6E

clearvars -except saveFigures
clc

load('./Data/Simulations/IndDepAnalysis/S_outparInd3_2_1')
% load('/Volumes/MELANOMA/Data/Data100')
load('./Data/Data100')
thres = Data100(:,1)./Data100(:,2).*Data100(:,5)*0.8;
param = 88;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
h = plot(855600:856600,S_outpar{1}(1:3,855600:856600),'LineWidth',2);
hold on 
plot(855600:856600,repmat(thres(param),1,length(855600:856600)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'}, {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255});
xlim([855600, 856600])
xticks([855600, 856600])
xticklabels({'0', '1000'})
yticks([0, 350])
ylim([0,350])
yticklabels({'0', '3.5'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
ylabel('gene product (10^2)')
xlabel('time')
box off

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS6/FigureS6E'])
    print('-dpdf',['./Figures/FigureS6/FigureS6E'])
end

%% SupplementeryFigure6F

clearvars -except saveFigures
clc

load('./Data/Simulations/3nodes/S_outpar3_2_1')
% load('/Volumes/MELANOMA/Data/Data100')
load('./Data/Data100')
thres = Data100(:,1)./Data100(:,2).*Data100(:,5)*0.8;
param = 88;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
h = plot(855800:856800,S_outpar{param}(1:3,855800:856800),'LineWidth',2);
hold on 
plot(855800:856800,repmat(thres(param),1,length(855800:856800)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'}, {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255});
xlim([855800, 856800])
xticks([855800 856800])
xticklabels({'0', '1000'})
yticks([0, 350])
ylim([0,350])
yticklabels({'0', '3.5'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
ylabel('gene product (10^2)')
xlabel('time')
box off

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS6/FigureS6F'])
    print('-dpdf',['./Figures/FigureS6/FigureS6F'])
end

%% SupplementeryFigure6G

clearvars -except saveFigures
clc

load('./Data/Simulations/IndDepAnalysis/S_outparDep3_2_1')
% load('/Volumes/MELANOMA/Data/Data100')
load('./Data/Data100')
thres = Data100(:,1)./Data100(:,2).*Data100(:,5)*0.8;
param = 88;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);
h = plot(855600:856600,S_outpar{1}(1:3,855600:856600),'LineWidth',2);
hold on 
plot(855600:856600,repmat(thres(param),1,length(855600:856600)),':','LineWidth',2, 'Color', 'k')
set(h, {'color'}, {[135 170 222]./255; [22 45 80]./255; [95 141 211]./255});
xlim([855600, 856600])
xticks([855600, 856600])
xticklabels({'0', '1000'})
yticks([0, 350])
ylim([0,350])
yticklabels({'0', '3.5'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
ylabel('gene product (10^2)')
xlabel('time')
box off

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS6/FigureS6G'])
    print('-dpdf',['./Figures/FigureS6/FigureS6G'])
end
