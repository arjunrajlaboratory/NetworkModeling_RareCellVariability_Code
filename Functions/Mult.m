%%
load('Data1000')
n_species = 5;
tcell = 1000;    %time per 'cell'
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load('S_outparSensAna26_3_2_1')

thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;

for isubnet = 100
    
    clear rare_par
    
    isubnet
    
    count = 1;
    
    loadsoutpar = sprintf('S_outparMult5_3_%d',isubnet);
    load(loadsoutpar)
    
    for i = 1:10
        
        param = 10*(isubnet-1)+i;
        threshold = thres(param,:);
        
        [maxjackpot_sol,desc_sol,rightskewed_sol,unimodal_sol,samp_sol,samp_time_sol,rand_time_sol] ...
            = analyzeQual_revision(n_species,tcell,threshold,S_outpar{i});
        
        sol{i}.maxjackpot = maxjackpot_sol;
        sol{i}.desc = desc_sol;
        sol{i}.rightskewed = rightskewed_sol;
        sol{i}.unimodal = unimodal_sol;
        sol{i}.samp = samp_sol;
        sol{i}.samp_time = samp_time_sol;
        sol{i}.time = rand_time_sol;
    end
    
    for j = 1:10
        if  sol{j}.maxjackpot == 1
            if sol{j}.desc == 1
                if sol{j}.rightskewed == 1
                    if sol{j}.unimodal == 1
                        rare_par(count) = 10*(isubnet-1)+j;
                        count = count + 1;
                    end
                end
            end
        end
    end
    
    save_sol = sprintf('sol1000Mult_%d_3_%d',n_species,isubnet);
    save(save_sol,'sol');
    
    if exist('rare_par') == 1
        save_rare = sprintf('rare_par1000Mult_%d_3_%d',n_species,isubnet);
        save(save_rare,'rare_par');
    else
        rare_par = zeros(0,1);
        save_rare = sprintf('rare_par1000Mult_%d_3_%d',n_species,isubnet);
        save(save_rare,'rare_par');
    end
end

%%
Rare = [];
for j = 1:100
    loadrare = sprintf('rare_par1000Mult_5_3_%d',j);
    load(loadrare)
    
    Rare = [Rare, rare_par];
end

%%
%%
load('sol1000Mult_5_3_97')
%number of highly expressed genes
figure
histogram(sol{8}.samp_time,'Facecolor',[0,0,0]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1);
set(gca,'linewidth',2)
box off
ylim([0,1])
yticks([0, 0.5, 1])
yticklabels({'0', '50', '100'})
xlim([-0.5,5.5])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
xlabel('no. highly expressed genes')
ylabel('% of cells')
title('Simulation')

%expression distribution
figure
subplot(2,1,1)
histogram(sol{8}.samp(1,:),'Facecolor',[135 170 222]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1, 'Binwidth',50);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([0, 0.5, 1])
yticklabels({'0', '50', '100'})
xticks([0, 350])
xlabel('gene product')
ylabel('% of cells')
title('Simulation')
xlim([0,350])
ylim([0 1])

%carpet
subplot(2,1,2)
plot(sol{8}.samp(1,:),ones(1000,1),'.','Color',[135 170 222]./255);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca, 'visible', 'off')
set(0, 'DefaultFigureRenderer', 'painters');
xlim([0,350])