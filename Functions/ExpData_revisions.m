%% Supplementary Figure (potentiall) 2

data = importfile('/Volumes/MELANOMA/Data/WM9_noDrug_20150618.txt');

for i = 12
% for i = [4,8,11,12,13,15,21,18,16]
    
    figure
    subplot(2,1,1)
    histogram(data(:,i),'Facecolor',[100,100,100]./255,'Normalization','probability',...
        'Binwidth', max(data(:,i))/10,'EdgeColor','none', 'FaceAlpha', 1);
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
    xlim([0,max(data(:,i))])
    
    %carpet
    subplot(2,1,2)
    plot(data(:,i),ones(length(data(:,12)),1),'.','Color',[100,100,100]./255);
    set(gca,'linewidth',2)
    box off
    set(0,'DefaultAxesFontName','Arial');
    set(gca,'FontSize',12)
    set(gca, 'visible', 'off')
    set(0, 'DefaultFigureRenderer', 'painters');
    xlim([0,max(data(:,i))])
    
    % if saveFigures
%     savepath = sprintf('/Volumes/MELANOMAII/Revisions/ExpData%d',i);
%     saveas(gcf, [savepath, '.fig']);
%     saveas(gcf, [savepath, '.pdf']);
%     saveas(gcf, [savepath, '.svg']);
    %     saveas(gcf,'./Figures/Figure2/Figure2C.fig');
    %     saveas(gcf,'./Figures/Figure2/Figure2C.pdf');
    %     saveas(gcf,'./Figures/Figure2/Figure2C.svg');
    % end
    
end

%% Heavy-tails for experimental data

data = importfile('/Volumes/MELANOMA/Data/WM9_noDrug_20150618.txt');

count = 1;

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

figure
plot(Com(:,3),Com(:,2),'.','MarkerSize',20, 'Color', [211,95,95]./255)
hold on
plot([1,300],[1,300],':','Color','k','LineWidth',2)
box off
set(gca,'linewidth',2)
xlabel('99th percentile of fitted exponential distribution')
ylabel({'99th percentile of expression'; 'distributions of simulations with'; 'rare coordinated high states'})
xticks([0,300])
yticks([0,300])
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
xlim([0,300])
ylim([0,300])

saveas(gcf,'/Volumes/MELANOMAII/Revisions/ExpDataHeavyTails.fig');
saveas(gcf,'/Volumes/MELANOMAII/Revisions/ExpDataHeavyTails.pdf');
saveas(gcf,'/Volumes/MELANOMAII/Revisions/ExpDataHeavyTails.svg');

%% heavy tails rest 

%get data from HeavytailsNew

figure
plot(Com(:,8),Com(:,7),'.','MarkerSize',20, 'Color', [211,95,95]./255)
hold on
plot([1,300],[1,300],':','Color','k','LineWidth',2)
box off
set(gca,'linewidth',2)
xlabel('99th percentile of fitted exponential distribution')
ylabel({'99th percentile of expression'; 'distributions of simulations with'; 'rare coordinated high states'})
xticks([0,300])
yticks([0,300])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
xlim([0,300])
ylim([0,300])

