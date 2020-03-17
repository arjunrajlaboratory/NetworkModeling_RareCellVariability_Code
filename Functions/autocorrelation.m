
for n_species = 3
    
    for inet = 2
        
        for isubnet = 1:10
            
            loadsol = sprintf('/Volumes/LEADE3/Data5nodes10000/Data%dnodes/S_outpar1000_%d_%d_%d',n_species,n_species,inet,isubnet);
            load(loadsol)
            
            for isoutpar = 1:length(S_outpar)
                [ta, a] = acf(S_outpar{isoutpar}(1,:)',1000-1);
                if isempty(find(ta<a,1)) == 0
                    T(countacf) = find(ta<a,1);
                else
                    T(countacf) = NaN;
                end
                countacf = countacf + 1;
            end
            
        end
    end
end
