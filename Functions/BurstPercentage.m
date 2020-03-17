Net = [2,4,10,80];
Subnet = [10,10,100,250];
Nspecies = [2,3,5,8];

count = 1;
countrare = 1;

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
                
                if n_species < 8
                    loadsol = sprintf('/Volumes/MELANOMA/Data/Simulations/%dnodes/S_outpar1000_%d_%d_%d.mat',n_species,n_species,inet,isubnet);
                    load(loadsol)
                else
                    if ismember(inet,[1,2,3,4,10:48]) == 1
                        loadsol = sprintf('/Volumes/MELANOMA/Data/Simulations/%dnodes/S_outpar1000%d_%d_%d.mat',n_species,n_species,inet,isubnet);
                        load(loadsol)
                    else
                        loadsol = sprintf('/Volumes/MELANOMAII/Data/Simulations/%dnodes/S_outpar1000%d_%d_%d.mat',n_species,n_species,inet,isubnet);
                        load(loadsol)
                    end
                end
                
                for param = rare_par_isubnet
                    
                    %number of species above threshold
                    f1 = sum(S_outpar{param-(isubnet-1)*length(S_outpar)}(1:n_species,1:end) > thres(param));
                    
                    %find all time points at which more than half of the species are above
                    %the threshold
                    if mod(n_species,2) == 1
                        L1 = find(f1 >= ceil(n_species/2));
                    else
                        L1 = find(f1 >= ceil(n_species/2) + 1);
                    end
                    
                    %find jackpot ends (difference greater than 1 between the time points
                    %showing jackpot)
                    diffJ = diff(L1);
                    TermJ_ind = find(diffJ > 1);
                    
                    L2 = zeros(1,0);
                    count_L = 1;
                    for iDiff = 1:length(TermJ_ind)
                        %if time period between the single jackpot regions is below 50,
                        %then take it as one jackpot
                        if diffJ(TermJ_ind(iDiff)) < 50
                            L2(count_L:count_L+diffJ(TermJ_ind(iDiff))-1) = L1(TermJ_ind(iDiff))+1:L1(TermJ_ind(iDiff)+1);
                            count_L = length(L2)+1;
                        end
                    end
                    
                    %L summarizes all time points at which the system is in a jackpot region
                    if isempty(L2) == 1
                        L = unique(L1);
                    else
                        L = unique([L1,L2]);
                    end
                    
                    %find jackpot times
                    L3 = find(f1 >= 1);
                    L3_start = L3([1,find(diff(L3)>50)+1]);
                    L3_end = L3([find(diff(L3)>50),length(L3)]);
                    
                    %find all time intervals in which at least one species is above threshold
                    %and which contains a jackpot in between start and end of jackpot region
                    countL = 1;
                    countLStart = 1;
                    for iL3 = 1:length(L3_start)
                        if isempty(intersect(L3_start(iL3):L3_end(iL3),L)) == 0
                            L(countL:countL+length(L3_start(iL3):L3_end(iL3))-1) =...
                                L3_start(iL3):L3_end(iL3);
                            LStart(countLStart) = L3_start(iL3);
                            LEnd(countLStart) = L3_end(iL3);
                            countLStart = countLStart + 1;
                            countL = length(L) + 1;
                        end
                    end
                    
                    TimeStartJ = LStart;
                    TimeFinishJ = LEnd;
                    
                    length_high(countrare) = length(unique(L));
                    bursts_high(countrare) = sum(S_outpar{param-(isubnet-1)*length(S_outpar)}(n_species+2,unique(L)));
                    
                    length_base(countrare) = length(S_outpar{param-(isubnet-1)*length(S_outpar)}(n_species+2,setdiff(1:1000000,unique(L))));
                    bursts_base(countrare) = sum(S_outpar{param-(isubnet-1)*length(S_outpar)}(n_species+2,setdiff(1:1000000,unique(L))));
                    
                    total_high(countrare) = sum(S_outpar{param-(isubnet-1)*length(S_outpar)}(n_species+2,1:end));
                    countrare = countrare + 1;
                end
                
            end
        end
    end
    count = count + 1;
end

[h,p,ks2stat] = kstest2((bursts_high./length_high),(bursts_base./length_base));

figure
[h,L,MX,MED]=violin([(bursts_high./length_high)',(bursts_base./length_base)'], 'facecolor', [150,150,150;150,150,150]./255);
xticks([1,2])
yticks([0,100])
xticklabels({'high time-region','baseline time-region'})
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
box off
set(gca,'linewidth',2)
ylabel('burst percentage')

