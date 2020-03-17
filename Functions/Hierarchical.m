load('/Volumes/MELANOMA/Data/Data1000')

S_outpar = generateDataHierarchical(1000,5,1000000,'no',Data1000);

%%
load('/Volumes/MELANOMA/Data/Data1000')
n_species = 5;
tcell = 1000;    %time per 'cell'
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('S_outparSensAna26_3_2_1')

thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;

for isubnet = 1:10
    
    clear rare_par
    
    isubnet
    
    count = 1;
    
    loadsoutpar = sprintf('/Volumes/MELANOMAII/Revisions/S_outparHierarchical5_1_%d',isubnet);
    load(loadsoutpar)
    
    for i = 1:100
        
        param = 100*(isubnet-1)+i;
        threshold = thres(param,:);
        
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
    
    for j = 1:100
        if  sol{j}.maxjackpot == 1
            if sol{j}.desc == 1
                if sol{j}.rightskewed == 1
                    if sol{j}.unimodal == 1
                        rare_par(count) = 100*(isubnet-1)+j;
                        count = count + 1;
                    end
                end
            end
        end
    end
    
    save_sol = sprintf('/Volumes/MELANOMAII/Revisions/sol1000Hierarchical_%d_1_%d',n_species,isubnet);
    save(save_sol,'sol');
    
    if exist('rare_par') == 1
        save_rare = sprintf('/Volumes/MELANOMAII/Revisions/rare_par1000Hierarchical_%d_1_%d',n_species,isubnet);
        save(save_rare,'rare_par');
    else
        rare_par = zeros(0,1);
        save_rare = sprintf('/Volumes/MELANOMAII/Revisions/rare_par1000Hierarchical_%d_1_%d',n_species,isubnet);
        save(save_rare,'rare_par');
    end
end

Rare = [];
for j = 1:10
    loadrare = sprintf('/Volumes/MELANOMAII/Revisions/rare_par1000Hierarchical_5_1_%d',j);
    load(loadrare)
    
    Rare = [Rare, rare_par];
end
