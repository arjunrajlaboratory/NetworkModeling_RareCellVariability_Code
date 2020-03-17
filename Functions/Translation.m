%%
load('Data1000')
n_species = 5;
tcell = 1000;    %time per 'cell'
loadsoutpar = sprintf('S_outparTransFaster5_3_1');
a = 10;
b = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('S_outparSensAna26_3_2_1')

load(loadsoutpar)
thres = Data1000(968,1)^2/Data1000(968,2)^2*Data1000(968,8)*0.8*a/b;
 
i = 1;
param = 1;
threshold = thres;

[maxjackpot_sol,desc_sol,rightskewed_sol,unimodal_sol,samp_sol,samp_time_sol,rand_time_sol] ...
    = analyzeQualTrans_revision(n_species,tcell,threshold,S_outpar{i});

sol{i}.maxjackpot = maxjackpot_sol;
sol{i}.desc = desc_sol;
sol{i}.rightskewed = rightskewed_sol;
sol{i}.unimodal = unimodal_sol;
sol{i}.samp = samp_sol;
sol{i}.samp_time = samp_time_sol;
sol{i}.time = rand_time_sol;
    
if  sol{i}.maxjackpot == 1
    if sol{i}.desc == 1
        if sol{i}.rightskewed == 1
            if sol{i}.unimodal == 1
                rare_par = i;
            end
        end
    end
end
   

%%
%number of highly expressed genes
figure
histogram(sol{1}.samp_time,'Facecolor',[0,0,0]./255,'Normalization','probability',...
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
histogram(sol{1}.samp(1,:),'Facecolor',[135 170 222]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1, 'Binwidth',50);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
yticks([0, 0.5, 1])
yticklabels({'0', '50', '100'})
xticks([0, 1400])
xlabel('gene product')
ylabel('% of cells')
title('Simulation')
xlim([0,1400])
ylim([0 1])

%carpet
subplot(2,1,2)
plot(sol{1}.samp(1,:),ones(1000,1),'.','Color',[135 170 222]./255);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca, 'visible', 'off')
set(0, 'DefaultFigureRenderer', 'painters');
xlim([0,1400])
