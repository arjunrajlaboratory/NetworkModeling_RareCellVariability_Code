
count = 1;
% for n_species = [2,3,5,8]
for n_species = 5
    
%     loadMiso = sprintf('/Volumes/MELANOMA/Data/M_iso%d',n_species);
%     load(loadMiso);
    load('M_iso5')
    
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
        
%         loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',n_species,n_species,inet);
%         load(loadrare)
%         
%         rare{count}(inet) = length(rare_par); 
    end
    count = count + 1;
end

Col{1} = [233,175,175]./255;
Col{2} = [211,95,95]./255;
Col{3} = [160,44,44]./255;
Col{4} = [120,33,33]./255;

figure
for iplot = 4:-1:1
    plot(L{iplot},rare{iplot}, '.','MarkerSize',20, 'Color', Col{iplot})
    hold on
end
box off
set(gca,'linewidth',2)
xlabel('characteristic distance')
ylabel({'number of simulations with rare coordinated high states'})
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
xlim([0,1])
ylim([0,70])
