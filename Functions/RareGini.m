Net = [2,4,10,80];
Subnet = [10,10,100,250];
Nspecies = [2,3,5,8];

countgini = 1;
countginiII = 1;
count = 1;

load('/Volumes/MELANOMA/Data/Data1000')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;

for n_species = [2,3,5,8]
    % for n_species = 8
    n_species
    for inet = 1:Net(count)
        inet
        
        loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',n_species,n_species,inet);
        load(loadrare)
        
        for isubnet = 1:Subnet(count)
            isubnet
            r1 = rare_par(rare_par > isubnet*1000/Subnet(count)-1000/Subnet(count));
            r2 = rare_par(rare_par <= isubnet*1000/Subnet(count));
            rare_par_isubnet = intersect(r1,r2);
            
            if isempty(rare_par_isubnet) == 0
                
                loadsol = sprintf('/Volumes/MELANOMAII/Data/CriteriaAnalysis/%dnodes/sol1000%d_%d_%d.mat',n_species,n_species,inet,isubnet);
                load(loadsol)
                
                for param = rare_par_isubnet
                    
                    genecount = sol{param-(isubnet-1)*1000/Subnet(count)}.samp(1,:);
                    
                    G(countgini) = gini(ones(1000,1),genecount);
                    countgini = countgini + 1;
                end
            end
            
            if isempty(rare_par_isubnet) == 1
                
                loadsol = sprintf('/Volumes/MELANOMAII/Data/CriteriaAnalysis/%dnodes/sol1000%d_%d_%d.mat',n_species,n_species,inet,isubnet);
                load(loadsol)
                
            end
            
            for paramII = setdiff(1000/Subnet(count)*(isubnet-1)+1:1000/Subnet(count)*isubnet,rare_par_isubnet)
                
                genecount = sol{paramII-(isubnet-1)*1000/Subnet(count)}.samp(1,:);
                
                GII(countginiII) = gini(ones(1000,1),genecount);
                countginiII = countginiII + 1;
            end
            
            
        end
    end
    count = count + 1;
end

figure
histogram(G,'Binwidth',0.025,'Normalization','probability','Facecolor',[100,100,100]./255,'Edgecolor','none')
hold on
histogram(GII,'Binwidth',0.025,'Normalization','probability','Facecolor',[200,200,200]./255,'Edgecolor','none')
alpha=1;
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