
count = 1;
for n_species = [2,3,5,8]
    
    loadmat = sprintf('/Volumes/MELANOMA/Data/M_iso%d',n_species);
    load(loadmat);
    
    for i = 1:length(M_iso)
        
        conn{count}(i) = sum(M_iso{i}(1,:));
        
        loops{count}(i) =  M_iso{i}(1,1);
        
        loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',n_species,n_species,i);
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
    
    count = count + 1;
end

%plot of network size
countplot = 1;
figure
for k = [2,3,5,8]
    plot(repmat(countplot,1,length(Ratio{countplot}(Ratio{countplot}>0))),...
        Ratio{countplot}(Ratio{countplot}>0),'.','MarkerSize',20, 'Color', 'k')
    hold on
    countplot = countplot + 1;
end
box off
set(gca,'linewidth',2)
xlabel('network size')
ylabel({'ratio of number of simulations with rare coordinated high states for networks without loops to networks with loops'})
xticks(1:4)
xticklabels({'2','3','5','8'})
yticks([0:5])
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
xlim([0,5])
ylim([0,5])

%plot of connectivity - network of size 8 only
figure
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
ylabel({'number of simulations with rare coordinated high states'})
xticks(2:2:16)
xticklabels({'1','2','3','4','5','6','7','8'})
yticks([0,35])
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
xlim([0,18])
ylim([0,35])
