countacf = 1;
countacf_empty = 1;

for n_species = 3

    for inet = 2
        
        loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',n_species,n_species,inet);
        load(loadrare)
        
        for isubnet = 9
            isubnet
            r1 = rare_par(rare_par > isubnet*1000/10-1000/10);
            r2 = rare_par(rare_par <= isubnet*1000/10);
            rare_par_isubnet = intersect(r1,r2);
            
            if isempty(rare_par_isubnet) == 0
                
                loadsol = sprintf('/Volumes/MELANOMA/Data/Simulations/%dnodes/S_outpar1000_%d_%d_%d.mat',n_species,n_species,inet,isubnet);
                load(loadsol)
                
                for param = rare_par_isubnet
                    [ta, a] = acf(S_outpar{param-(isubnet-1)*length(S_outpar)}(1,:)',1000-1);
                    if isempty(find(ta<a,1)) == 0
                        T(countacf) = find(ta<a,1);
                    else
                        T(countacf) = NaN;
                        End(countacf_empty) = ta(end);
                        countacf_empty = countacf_empty + 1;
                    end
                    countacf = countacf + 1;
                end
                
            end
        end
    end
end
