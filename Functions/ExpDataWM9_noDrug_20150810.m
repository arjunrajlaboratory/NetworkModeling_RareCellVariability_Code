data = importfile('/Users/lea.schuh/Documents/PhD/ICB/Melanoma/Revisions/WM9_noDrug_20150810.txt');

for i = 12
    
    figure
    subplot(2,1,1)
    histogram(data(:,i),'Facecolor',[100,100,100]./255,'Normalization','probability',...
        'Binwidth', 50,'EdgeColor','none', 'FaceAlpha', 1);
    set(gca,'linewidth',2)
    box off
    set(0,'DefaultAxesFontName','Arial');
    set(gca,'FontSize',12)
    yticks([0, 0.5, 1])
    yticklabels({'0', '50', '100'})
    xlabel('gene product')
    ylabel('% of cells')
    title('Data')
    ylim([0 1])
    xlim([0,750])
    
    %carpet
    subplot(2,1,2)
    plot(data(:,i),ones(length(data(:,12)),1),'.','Color',[100,100,100]./255);
    set(gca,'linewidth',2)
    box off
    set(0,'DefaultAxesFontName','Arial');
    set(gca,'FontSize',12)
    set(gca, 'visible', 'off')
    set(0, 'DefaultFigureRenderer', 'painters');
    xlim([0,750])
    
end