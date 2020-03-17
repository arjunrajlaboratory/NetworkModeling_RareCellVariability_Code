
addpath(genpath(pwd));

saveFigures = true;

%% SupplementeryFigure9A

clearvars -except saveFigures
clc

% load('/Volumes/MELANOMA/Data/M_iso2')
load('./Data/M_iso2')
count = 0;

for j = 1:length(M_iso)
    
    pos = [0.1,0.32,0.54,0.76];
    xlen = 0.2;
    ylen = 0.2;
    
    figure
    i = 1;
    subplot('Position',[pos(i),pos(1),xlen,ylen]);
    
    count = count + 1;
    G = digraph(M_iso{j});
    h = plot(G,'MarkerSize',5, 'NodeLabel',{}, 'LineWidth', 1, 'ArrowSize', 5, 'Edgecolor','k','Nodecolor',[150,150,150]./255,'EdgeAlpha',1);
    ax = gca;
    ax.Visible = 'off';
    
    if saveFigures
        set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%         path = sprintf('/Volumes/MELANOMAII/Figures/FigureS9/FigureS9A%d',j);
         path = sprintf('./Figures/FigureS9/FigureS9A%d',j);
        print('-dpdf',[path])
    end
end


%% SupplementeryFigure9B

clearvars -except saveFigures
clc

% load('/Volumes/MELANOMA/Data/M_iso3')
load('./Data/M_iso3')
count = 0;

for j = 1:length(M_iso)
    
    pos = [0.1,0.32,0.54,0.76];
    xlen = 0.2;
    ylen = 0.2;
    
    figure
    i = 1;
    subplot('Position',[pos(i),pos(1),xlen,ylen]);
    
    count = count + 1;
    G = digraph(M_iso{j});
    h = plot(G,'MarkerSize',5, 'NodeLabel',{}, 'LineWidth', 1, 'ArrowSize', 5, 'Edgecolor','k','Nodecolor',[150,150,150]./255,'EdgeAlpha',1);
    ax = gca;
    ax.Visible = 'off';
    
    if saveFigures
        set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%         path = sprintf('/Volumes/MELANOMAII/Figures/FigureS9/FigureS9B%d',j);
        path = sprintf('./Figures/FigureS9/FigureS9B%d',j);
        print('-dpdf',[path])
    end
end


%% SupplementeryFigure9C

clearvars -except saveFigures
clc

% load('/Volumes/MELANOMA/Data/M_iso5')
load('./Data/M_iso5')
count = 0;

for j = 1:length(M_iso)
    
    pos = [0.1,0.32,0.54,0.76];
    xlen = 0.2;
    ylen = 0.2;
    
    figure
    i = 1;
    subplot('Position',[pos(i),pos(1),xlen,ylen]);
    
    count = count + 1;
    G = digraph(M_iso{j});
    h = plot(G,'MarkerSize',5, 'NodeLabel',{}, 'LineWidth', 1, 'ArrowSize', 5, 'Edgecolor','k','Nodecolor',[150,150,150]./255,'EdgeAlpha',1);
    ax = gca;
    ax.Visible = 'off';
    
    if saveFigures
        set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%         path = sprintf('/Volumes/MELANOMAII/Figures/FigureS9/FigureS9C%d',j);
        path = sprintf('./Figures/FigureS9/FigureS9C%d',j);
        print('-dpdf',[path])
   end
    
end

%% SupplementeryFigure9D

clearvars -except saveFigures
clc

% load('/Volumes/MELANOMA/Data/M_iso8')
load('./Data/M_iso8')
count = 0;

for j = 1:length(M_iso)
    
    pos = [0.1,0.32,0.54,0.76];
    xlen = 0.2;
    ylen = 0.2;
    
    figure
    i = 1;
    subplot('Position',[pos(i),pos(1),xlen,ylen]);
    
    count = count + 1;
    G = digraph(M_iso{j});
    
    h = plot(G,'MarkerSize',5, 'NodeLabel',{}, 'LineWidth', 1, 'ArrowSize', 5, 'Edgecolor','k','Nodecolor',[150,150,150]./255,'EdgeAlpha',1);
    ax = gca;
    ax.Visible = 'off';
    
    if saveFigures
        set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%         path = sprintf('/Volumes/MELANOMAII/Figures/FigureS9/FigureS9D%d',j);
        path = sprintf('./Figures/FigureS9/FigureS9D%d',j);
        print('-dpdf',[path])
    end
    
end