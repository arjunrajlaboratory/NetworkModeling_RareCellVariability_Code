countrare = 1;

load('/Volumes/MELANOMA/Data/Data1000')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;

% for n_species = [2,3,5,8]
for n_species = 3
    
    for inet = 2
        
        loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',n_species,n_species,inet);
        load(loadrare)
        
        for isubnet = 1:10
            
            isubnet
            
            r1 = rare_par(rare_par > isubnet*1000/10-1000/10);
            r2 = rare_par(rare_par <= isubnet*1000/10);
            rare_par_isubnet = intersect(r1,r2);
            
            if isempty(rare_par_isubnet) == 0
                
                %                 loadsol = sprintf('/Volumes/LEADE1/Paper_Code/S_outparInd1000%d_%d_%d.mat',n_species,inet,isubnet);
                %                 load(loadsol)
                loadsol = sprintf('/Volumes/MELANOMA/Data/Simulations/%dnodes/S_outpar1000_%d_%d_%d.mat',n_species,n_species,inet,isubnet);
                load(loadsol)
                
                for param = rare_par_isubnet
                    
                    clearvars -except countrare Data1000 thres n_species inet rare_par isubnet rare_par_isubnet length_L_32 S_outpar param length_L_32_dep
                    
                    %number of species above threshold
                    f1 = S_outpar{param-(isubnet-1)*length(S_outpar)}(1,1:end) > thres(param);
                 
                    %find jackpot times
                    L3 = find(f1 == 1);
                    if isempty(L3) == 0
                        L3_start = L3([1,find(diff(L3)>1)+1]);
                        L3_end = L3([find(diff(L3)>1),length(L3)]);

                        for i = 1:length(L3_start)
                            length_L_32_dep{countrare}(i) = length(L3_start(i):L3_end(i));
                        end
                        countrare = countrare + 1;
                    else
                        length_L_32_dep{countrare} = NaN;
                        countrare = countrare + 1;
                    end
                end
                
            end
        end
    end
end

%% figure
load('length_L_32_dep')
load('length_L_32')

for i = 1:26
    Sdep(i) = sum(length_L_32_dep{i});
    S(i) = sum(length_L_32{i});
    Ndep(i) = length(length_L_32_dep{i});
    N(i) = length(length_L_32{i});
end


figure
[h,L,MX,MED]=violin((Ndep./N)', 'facecolor', [150,150,150;150,150,150]./255);
hold on
plot(ones(length((Ndep./N)),1)+1/50*randn(length((Ndep./N)),1),(Ndep./N),'.','Markersize', 20, 'Color', 'k')
hold on 
plot([0,2],[1,1],'-')
xticks([1,2])
% yticks([0,100])
xticklabels({'network','indepedent'})
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
box off
set(gca,'linewidth',2)
ylabel('number of high states')


figure
[h,L,MX,MED]=violin((Sdep./S)', 'facecolor', [150,150,150;150,150,150]./255);
hold on
plot(ones(length((Sdep./S)),1)+1/50*randn(length((Sdep./S)),1),(Sdep./S),'.','Markersize', 20, 'Color', 'k')
hold on 
plot([0,2],[1,1],'-')
xticks([1,2])
% yticks([0,100])
xticklabels({'network','indepedent'})
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
box off
set(gca,'linewidth',2)
ylabel('number of high states')

figure
plot(ones(1,26),S./Sdep, '.','Markersize', 20, 'Color', 'k');
hold on
plot(repmat(2,1,26),N./Ndep,'.', 'Markersize', 20, 'Color', 'k');
xticks([1,2])
yticks([0,0.5,1])
xticklabels({'total time in high state','number of high states'})
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
box off
set(gca,'linewidth',2)
ylabel('ratio of independent to network')
xlim([0,3])
ylim([0,1.2])

figure
[h,L,MX,MED]=violin([Ndep',N'], 'facecolor', [150,150,150;150,150,150]./255);
xticks([1,2])
% yticks([0,100])
xticklabels({'network','indepedent'})
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
box off
set(gca,'linewidth',2)
ylabel('number of high states')


