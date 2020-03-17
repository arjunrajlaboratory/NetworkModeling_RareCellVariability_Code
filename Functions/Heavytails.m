n_species = 3;
Net = 1;
Subnet = 10;


count = 1;

for inet = 1:Net
    
    loadrare = sprintf('/Volumes/LEADE1/Paper_Code/%dnodes/rare_par10003_%d_%d.mat',n_species,n_species,inet);
%     loadrare = sprintf('/Volumes/LEADE3/Data5nodes10000/rare_par10003_%d_%d.mat',n_species,inet);
    load(loadrare)
    
    for isubnet = 1:Subnet
        
        r1 = rare_par(rare_par > isubnet*1000/Subnet-1000/Subnet);
        r2 = rare_par(rare_par <= isubnet*1000/Subnet);
        rare_par_isubnet = intersect(r1,r2);
        
        if isempty(rare_par_isubnet) == 0
            
%             loadsol = sprintf('/Volumes/LEADE1/Paper_Code/%dnodes/sol10003%d_%d_%d.mat',n_species,n_species,inet,isubnet);
            loadsol = sprintf('/Volumes/MELANOMA/Data/CriteriaAnalysis/%dnodes/Replicates3/sol10003%d_%d_%d.mat',n_species,n_species,inet,isubnet);
            load(loadsol)
            
            for i = rare_par_isubnet
                DataSim = sol{i-(length(sol)*isubnet-length(sol))}.samp(1,:);
                PrctlSim = prctile(sol{i-(length(sol)*isubnet-length(sol))}.samp(1,:),99);

                a = histcounts(DataSim,'BinWidth',1);
                bin = 1;
                
%                 if find(a == max(a)) <= 15

                    anew = max(a);
                    error = 100;
                    countwhile = 0;
                    while abs(error) > 1
                        if error > 0
                            anew = anew + 10;
                        else
                            anew = anew - 10;
                            if anew < 0
                                DataExp = zeros(0,1);
                                break
                            end
                        end
                        DataExp = exprnd(anew,1,1000);
                        b = histcounts(DataExp,'BinWidth',bin);
                        error = max(b)-max(a);
                        countwhile = countwhile + 1;
                        if countwhile > 1000
                            DataExp = zeros(0,1);
                            break
                        end
                    end

                    A(count) = anew;
                    if isempty(DataExp) == 0
                        
                        PrctlExp = prctile(DataExp,99); 
                        Com(count,:) = [1,n_species,inet,isubnet,i-(length(sol)*isubnet-length(sol)),i,PrctlSim,PrctlExp];
                        count = count+1;

                    else
                        Com(count,:) = [1,n_species,inet,isubnet,i-(length(sol)*isubnet-length(sol)),i,NaN,NaN]; %0 for not run through
                        count = count+1;
                    end
                    
%                 end
                
                
            end
        end
    end
end
