
addpath(genpath(pwd));

saveFigures = true;

%% SupplementeryFigure5A

%normlaizing plot Fig 2D
Net = [2,4,10,80];
ncount = 0;

for n_species = [2,3,5,8]
    ncount = ncount + 1;
    
    R1 = zeros(1,0);
    R2 = zeros(1,0);
    R3 = zeros(1,0);
   
%     loadMiso = sprintf('/Volumes/MELANOMA/Data/M_iso%d',n_species);
    loadMiso = sprintf('./Data/M_iso%d',n_species);
    load(loadMiso)
    
    for inet = 1:Net(ncount)
        
        
        rare_par = [];
        
%         loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',...
%             n_species,n_species,inet);
        loadrare = sprintf('./Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',...
            n_species,n_species,inet);
        load(loadrare)
        R1 = [R1,rare_par];
        
        
        rare_par = [];
        
%         loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/Replicates2/rare_par10002_%d_%d.mat',...
%             n_species,n_species,inet);
        loadrare = sprintf('./Data/RareParameters/%dnodes/Replicates2/rare_par10002_%d_%d.mat',...
            n_species,n_species,inet);
        load(loadrare)
        R2 = [R2,rare_par];
        
        
        rare_par = [];
        
%         loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/Replicates3/rare_par10003_%d_%d.mat',...
%             n_species,n_species,inet);
        loadrare = sprintf('./Data/RareParameters/%dnodes/Replicates3/rare_par10003_%d_%d.mat',...
            n_species,n_species,inet);
        load(loadrare)
        R3 = [R3,rare_par];
        
    end
    
    RAll(ncount,:) = [length(R1),length(R2),length(R3)];
end

%normalized by network size
pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);

boxplot([RAll(1,:)'./2,RAll(2,:)'./3, RAll(3,:)'./5, RAll(4,:)'./8],'Colors', [0,0,0]./255)
ylim([0,60])
xlabel('network size')
xticklabels({'2','3','5','8'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
ylabel('simulations')
yticks([0,30,60])
box off
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
set(gca,'linewidth',2)

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS5/FigureS5A'])
    print('-dpdf',['./Figures/FigureS5/FigureS5A'])
end

%% SupplementeryFigure5B

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);

boxplot([RAll(1,:)'./2,RAll(2,:)'./4, RAll(3,:)'./10, RAll(4,:)'./80],'Colors', [0,0,0]./255)
ylim([0,60])
xlabel('network size')
xticklabels({'2','3','5','8'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
ylabel('simulations')
yticks([0,30,60])
box off
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
set(gca,'linewidth',2)

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS5/FigureS5B'])
    print('-dpdf',['./Figures/FigureS5/FigureS5B'])
end

%% SupplementeryFigure5C

clearvars -except saveFigures
clc

count = 0;
Net = [2,4,80];

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure

for n_species = [2,3,8]
    count = count +1;
    for inet = 1:Net(count)
        rare_par = [];
%         loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',...
%             n_species,n_species,inet);
        loadrare = sprintf('./Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',...
            n_species,n_species,inet);
        load(loadrare)
        R(inet) = length(rare_par);
    end

%     loadmat = sprintf('/Volumes/MELANOMA/Data/M_iso%d',n_species);
    loadmat = sprintf('./Data/M_iso%d',n_species);
    load(loadmat)
    for jnet = 1:length(M_iso)
        C(jnet) = sum(M_iso{jnet}(:,1));
    end

    j = count;
    subplot('Position',[pos(j),pos(1),xlen,ylen]);
    plot(C,R,'.','Markersize',20,'Color', [0, 0, 0]./255)
    xticks([1:n_species])
    ylim([0,70])
    xlim([0.5,n_species+0.5])
    box off
    set(0,'DefaultAxesFontName','Arial'); 
    set(gca,'FontSize',12)
    set(gca,'linewidth',2)
    if count == 1
        ylabel({'simulations with rare';'coordinated high states'})
        yticks([0,35,70])
    else
        yticks([])
    end
    xlabel('connectivity')
end

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS5/FigureS5C'])
    print('-dpdf',['./Figures/FigureS5/FigureS5C'])
end

%% SupplementeryFigure5D

clearvars -except saveFigures
clc

count = 0;
Net = [2,4,10,80];
Col{1} = [233,175,175]./255;
Col{2} = [211,95,95]./255;
Col{3} = [160,44,44]./255;
Col{4} = [120,33,33]./255;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);

for n_species = [2,3,5,8]
    clearvars Cjitter
    count = count +1;
    for inet = 1:Net(count)
        rare_par = [];
%         loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',...
%             n_species,n_species,inet);
        loadrare = sprintf('./Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',...
            n_species,n_species,inet);
        load(loadrare)
        R(inet) = length(rare_par);
    end

%     loadmat = sprintf('/Volumes/MELANOMA/Data/M_iso%d',n_species);
    loadmat = sprintf('./Data/M_iso%d',n_species);
    load(loadmat)
    for jnet = 1:length(M_iso)
        C(jnet) = sum(M_iso{jnet}(:,1));
    end
    
    Cjitter = C+1/10*randn(length(C),1)';
    
    plot(Cjitter,R,'.','Markersize',10,'Color', Col{count})
    hold on
end

xticks(1:8)
yticks([0,35,70])
ylim([0,70])
xlim([0.5,8.5])
box off
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
set(gca,'linewidth',2)
xlabel('connectivity')
ylabel({'simulations with rare';'coordinated high states'})

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS5/FigureS5D'])
    print('-dpdf',['./Figures/FigureS5/FigureS5D'])
end


%% SupplementeryFigure5E

count = 1;
for n_species = [2,3,5,8]
    
%     loadmat = sprintf('/Volumes/MELANOMA/Data/M_iso%d',n_species);
    loadmat = sprintf('./Data/M_iso%d',n_species);
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
    
    Ratio{count} = rare{count}(loops{count} == 0)./rare{count}(loops{count} == 1);
    
    if n_species == 8
        for j = unique(conn{count})
            indconn = find(conn{count} == j);
            indnoloops = find(loops{count} == 0);
            indloops = find(loops{count} == 1);
            
            indconnnoloops = intersect(indconn,indnoloops);
            indconnloops = intersect(indconn,indloops);
            
            for m = 1:length(indconnnoloops)
                Rare_noloops{j}(m) = length(Rare_par{indconnnoloops(m)});
            end
            for n = 1:length(indconnloops)
                Rare_loops{j}(n) = length(Rare_par{indconnloops(n)});
            end
        end
    end
    
    if n_species == 5
        for j = unique(conn{count})
            indconn5 = find(conn{count} == j);
            indnoloops5 = find(loops{count} == 0);
            indloops5 = find(loops{count} == 1);
            
            indconnnoloops5 = intersect(indconn5,indnoloops5);
            indconnloops5 = intersect(indconn5,indloops5);
            
            for m = 1:length(indconnnoloops5)
                Rare_noloops5{j}(m) = length(Rare_par{indconnnoloops5(m)});
            end
            for n = 1:length(indconnloops5)
                Rare_loops5{j}(n) = length(Rare_par{indconnloops5(n)});
            end
        end
    end
    
    count = count + 1;
end

%plot of connectivity - network of size 8 only
pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);

for h = unique(conn{3})
   plot(repmat(2*h-0.25,1,length(Rare_loops5{h})), Rare_loops5{h}, '.','MarkerSize',20, 'Color', 'k')
   hold on
   if h < 5
        plot(repmat(2*h+0.25,1,length(Rare_noloops5{h})), Rare_noloops5{h}, '.','MarkerSize',20, 'Color', [100,100,100]./255)
   end
end
box off
set(gca,'linewidth',2)
xlabel('connectivity')
ylabel({'simulations with rare';'coordinated high states'})
xticks(2:2:10)
xticklabels({'1','2','3','4','5'})
yticks([0,35])
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
xlim([0,12])
ylim([0,60])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS5/FigureS5Eleft'])
    print('-dpdf',['./Figures/FigureS5/FigureS5Eleft'])
end

%plot of connectivity - network of size 8 only
pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);

for h = unique(conn{4})
   plot(repmat(2*h-0.25,1,length(Rare_loops{h})), Rare_loops{h}, '.','MarkerSize',20, 'Color', 'k')
   hold on
   if h < 8
        plot(repmat(2*h+0.25,1,length(Rare_noloops{h})), Rare_noloops{h}, '.','MarkerSize',20, 'Color', [100,100,100]./255)
   end
end
box off
set(gca,'linewidth',2)
xlabel('connectivity')
ylabel({'simulations with rare';'coordinated high states'})
xticks(2:2:16)
xticklabels({'1','2','3','4','5','6','7','8'})
yticks([0,35])
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
xlim([0,18])
ylim([0,35])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS5/FigureS5Eright'])
    print('-dpdf',['./Figures/FigureS5/FigureS5Eright'])
end

%% SupplementeryFigure5F

count = 1;
for n_species = [2,3,5,8]
    
%     loadMiso = sprintf('/Volumes/MELANOMA/Data/M_iso%d',n_species);
    loadMiso = sprintf('./Data/M_iso%d',n_species);
    load(loadMiso);
    
    for inet = 1:length(M_iso)
    
        n = size(M_iso{inet},1);
        
        for in = 1:n
            for jn = 1:n
                D{count}{inet}(in,jn) = length(shortestpath(digraph(M_iso{inet}),in,jn))-1;
            end
        end
        
        for in= 1:n
            if M_iso{inet}(in,in) == 1
                D{count}{inet}(in,in) = D{count}{inet}(in,in)+1;
            else
                ind = find(M_iso{inet}(:,in)==1);
                I = [];
                for iind = ind'
                    I = [I,length(shortestpath(digraph(M_iso{inet}),in,iind))-1];
                end
                D{count}{inet}(in,in) = D{count}{inet}(in,in) + min(I) + 1;
            end
        end

        L{count}(inet) = mean(mean(D{count}{inet}))/n;
        
%         loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',...
%             n_species,n_species,inet);
        loadrare = sprintf('./Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',...
            n_species,n_species,inet);
        load(loadrare)
        
        rare{count}(inet) = length(rare_par); 
    end
    count = count + 1;
end

Col{1} = [233,175,175]./255;
Col{2} = [211,95,95]./255;
Col{3} = [160,44,44]./255;
Col{4} = [120,33,33]./255;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

figure
j = 1;
subplot('Position',[pos(j),pos(1),xlen,ylen]);

for iplot = 4:-1:1
    plot(L{iplot},rare{iplot}, '.','MarkerSize',10, 'Color', Col{iplot})
    hold on
end
box off
set(gca,'linewidth',2)
xlabel('characteristic distance')
ylabel({'simulations'})
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
xlim([0,1])
ylim([0,70])
yticks([0,35,70])
xticks([0:0.2:1])

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
%     print('-dpdf',['/Volumes/MELANOMAII/Figures/FigureS5/FigureS5F'])
    print('-dpdf',['./Figures/FigureS5/FigureS5F'])
end
