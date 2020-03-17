%find out all stable low/stable high parameters
% load('/Volumes/MELANOMAII/Data/classAll2.mat')
load('classAll2.mat')
C2_low = class(find(class(:,4) == 0),3);
C2_high = class(find(class(:,4) == 3),3);

% load('/Volumes/MELANOMAII/Data/classAll3.mat')
load('classAll3.mat')
C3_low = class(find(class(:,4) == 0),3);
C3_high = class(find(class(:,4) == 3),3);

% load('/Volumes/MELANOMAII/Data/classAll5.mat')
load('classAll5.mat')
C5_low = class(find(class(:,4) == 0),3);
C5_high = class(find(class(:,4) == 3),3);

% load('/Volumes/MELANOMAII/Data/classAll8.mat')
load('classAll8.mat')
C8_low = class(find(class(:,4) == 0),3);
C8_high = class(find(class(:,4) == 3),3);

CAll_low = [C2_low;C3_low;C5_low;C8_low];
CAll_high = [C2_high;C3_high;C5_high;C8_high];

for i = 1:1000
    H_low{i} = find(CAll_low == i);
    L_low(i) = length(find(CAll_low == i));
    
    H_high{i} = find(CAll_high == i);
    L_high(i) = length(find(CAll_high == i));
end

ParLow = find(L_low == max(L_low));
ParHigh = find(L_high == max(L_high));

%8 jackpot parameter sets, 2 stable low, 2 stable high
PAll = [26 92 133 183 544 702 915 968 504 889 125 944];

% load('/Volumes/MELANOMA/Data/Data1000.mat')
load('Data1000')

Data12 = Data1000(PAll,:);

%for networks connectivity 1, 6, 8, 10
% load('/Volumes/MELANOMA/Data/M_iso10.mat')
load('M_iso10')
for j = 1:length(M_iso)
    Con(j) = sum(M_iso{j}(1,:));
end

Con1 = find(Con == 1);
Con10 = find(Con == 10);
Con6 = find(Con == 6);
Con8 = find(Con == 8);

%networks: connectivity 1, 10 and two of each connectivity, 6 and 8
NetworkAll = [Con1,Con6([1,end]),Con8([1,end]),Con10];

generateData10(12,10,1000000,'no',Data12,NetworkAll)

%%
%load data sets (run on server)

load('/Volumes/MELANOMA/Data/Data1000.mat')
PAll = [26 92 133 183 544 702 915 968 504 889 125 944];

Data12 = Data1000(PAll,:);

load('/Volumes/MELANOMA/Data/M_iso10.mat')
for j = 1:length(M_iso)
    Con(j) = sum(M_iso{j}(1,:));
end

Con1 = find(Con == 1);
Con10 = find(Con == 10);
Con6 = find(Con == 6);
Con8 = find(Con == 8);

%networks: connectivity 1, 10 and two of each connectivity, 6 and 8
NetworkAll = [Con1,Con6([1,end]),Con8([1,end]),Con10];

for inet = NetworkAll
    
    loadsol = sprintf('/Volumes/MELANOMAII/Revisions/S_outpar1210_%d_1',inet);
    load(loadsol)
    
    n_species = 10;
    tcell = 1000;    %time per 'cell'
    doplot = 'no';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % load('S_outparSensAna26_3_2_1')
    count = 1;
    thres = Data12(:,1)./Data12(:,2).*Data12(:,8)*0.8;
    
    
    for i = 1:length(PAll)
        
        clear rare_par
        
        param = i;
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
    
    for j = 1:length(PAll)
        if  sol{j}.maxjackpot == 1
            if sol{j}.desc == 1
                if sol{j}.rightskewed == 1
                    if sol{j}.unimodal == 1
                        rare_par(count) = j;
                        count = count + 1;
                    end
                end
            end
        end
    end
    
    save_sol = sprintf('/Volumes/MELANOMAII/Revisions/sol1000_%d_%d',n_species,inet);
    save(save_sol,'sol');
    
    if exist('rare_par') == 1
        save_rare = sprintf('/Volumes/MELANOMAII/Revisions/rare_par1000_%d_%d',n_species,inet);
        save(save_rare,'rare_par');
    else
        rare_par = zeros(0,1);
        save_rare = sprintf('/Volumes/MELANOMAII/Revisions/rare_par1000_%d_%d',n_species,inet);
        save(save_rare,'rare_par');
    end
    
end

%%

load('/Volumes/MELANOMA/Data/M_iso10.mat')
for j = 1:length(M_iso)
    Con(j) = sum(M_iso{j}(1,:));
end

Con1 = find(Con == 1);
Con10 = find(Con == 10);
Con6 = find(Con == 6);
Con8 = find(Con == 8);

%networks: connectivity 1, 10 and two of each connectivity, 6 and 8
NetworkAll = [Con1,Con6([1,end]),Con8([1,end]),Con10];

count = 1;

for inet = NetworkAll
    
    loadsol = sprintf('/Volumes/MELANOMAII/Revisions/rare_par1000_10_%d',inet);
    load(loadsol)
    
    Rare{count} = rare_par;
    count = count + 1;
    
    clear rare_par
    
end

%% plot

loadsoutpar = sprintf('/Volumes/MELANOMAII/Revisions/S_outpar1210_1_1');
load(loadsoutpar)

loadsol = sprintf('/Volumes/MELANOMAII/Revisions/sol1000_10_1');
load(loadsol)
    
load('/Volumes/MELANOMA/Data/Data1000.mat')
PAll = [26 92 133 183 544 702 915 968 504 889 125 944];

Data12 = Data1000(PAll,:);
thres = Data12(:,1)./Data12(:,2).*Data12(:,8)*0.8;

%colors?
figure
plot(0:300,S_outpar{8}(1:10,447100:447400),'LineWidth',4)
hold on
plot([0,300],[thres(8),thres(8)],':','LineWidth',...
    4, 'Color', 'k')
yticks([0, 350])
yticklabels({'0', '3.5'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12,'linewidth',2)
ylabel('gene product (10^2)')
xlabel('time')
box off
legend

%right panel 'Simulation'
figure
param = 8;
histogram(sol{8}.samp_time,'Facecolor',[0,0,0]./255,'Normalization','probability',...
    'EdgeColor','none', 'FaceAlpha', 1);
set(gca,'linewidth',2)
box off
ylim([0,1])
yticks([0, 0.5, 1])
yticklabels({'0', '50', '100'})
xlim([-0.5,10.5])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
xlabel('no. highly expressed genes')
ylabel('% of cells')
title('Simulation')

%right panel 'Simulation'
figure
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
figure
plot(sol{8}.samp(1,:),ones(1000,1),'.','Color',[135 170 222]./255);
set(gca,'linewidth',2)
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca, 'visible', 'off')
set(0, 'DefaultFigureRenderer', 'painters');
xlim([0,350])

