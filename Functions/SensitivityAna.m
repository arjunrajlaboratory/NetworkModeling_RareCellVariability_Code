%sensitivity analysis
%for network architecture 3.2

%new parameter sets

clearvars -except saveFigures
clc

% JP = zeros(1,0);
%
% for n_species = [2,3,5,8]
%
%     loadMiso = sprintf('/Volumes/MELANOMA/Data/M_iso%d',n_species);
%     load(loadMiso);
%
%     for inet = 1:length(M_iso)
%         loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/rare_par1000_%d_%d',n_species,n_species,inet);
%         load(loadrare)
%         if isempty(rare_par) == 0
%             JP = [JP,rare_par];
%         end
%     end
% end
%
% for inum = 1:1000
%     ParFreq(inum) = length(find(JP == inum));
% end
%
% Par_Original = find(ParFreq == max(ParFreq));

load('/Volumes/MELANOMA/Data/Data1000')
ParSet_Original = Data1000(92,:);

%get min/max
% MinPar = min(Data1000([26 92 133 183 544 702 915 968],[3,5,6]));
% MaxPar = max(Data1000([26 92 133 183 544 702 915 968],[3,5,6]));

%build new parameter sets
DataSens = repmat(ParSet_Original,70,1);
DataSens(1:10,1) = linspace(0.01,1,10);
DataSens(11:20,2) = linspace(0.001,0.1,10);
DataSens(21:30,3) = linspace(0.001,0.1,10);
DataSens(31:40,4) = linspace(0.1,10,10);
DataSens(41:50,5) = linspace(0.1,1,10);
DataSens(51:60,6) = linspace(0.01,0.1,10);
DataSens(61:70,8) = linspace(2,100,10);


DataSens(:,7) = 0.95*DataSens(:,1)./DataSens(:,2).*DataSens(:,8);

S_outpar = generateDataSensAna(70,3,1000000,'no',DataSens);

%analyze the resulting simulations

n_species = 3;
inet = 2;
tcell = 1000;    %time per 'cell'
doplot = 'no';
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('S_outparSensAna26_3_2_1')

count = 1;

thres = DataSens(:,1)./DataSens(:,2).*DataSens(:,8)*0.8;


for i = 1:70
    
    param = i;
    threshold = thres(param,:);
    
    [maxjackpot_sol,desc_sol,rightskewed_sol,unimodal_sol,samp_sol,samp_time_sol,rand_time_sol] ...
        = analyzeQual(n_species,tcell,threshold,S_outpar{i},doplot,param,inet);
    
    sol{i}.maxjackpot = maxjackpot_sol;
    sol{i}.desc = desc_sol;
    sol{i}.rightskewed = rightskewed_sol;
    sol{i}.unimodal = unimodal_sol;
    sol{i}.samp = samp_sol;
    sol{i}.samp_time = samp_time_sol;
    sol{i}.time = rand_time_sol;
end

for j = 1:70
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

save_sol = sprintf('/Volumes/MELANOMAII/Revisions/sol1000SensAna92_3%d_%d',n_species,inet);
save(save_sol,'sol');

if exist('rare_par') == 1
    save_rare = sprintf('/Volumes/MELANOMAII/Revisions/rare_par1000SensAna92_3_%d_%d',n_species,inet);
    save(save_rare,'rare_par');
else
    rare_par = zeros(0,1);
    save_rare = sprintf('/Volumes/MELANOMAII/Revisions/rare_par1000SensAna92_3_%d_%d',n_species,inet);
    save(save_rare,'rare_par');
end

load('/Volumes/MELANOMA/Data/Data1000')
I = [26 92 133 183 544 702 915 968];
DataSens = repmat(Data1000(I(1),:),70,1);
DataSens(1:10,1) = linspace(0.01,1,10);
DataSens(11:20,2) = linspace(0.001,0.1,10);
DataSens(21:30,3) = linspace(0.001,0.1,10);
DataSens(31:40,4) = linspace(0.1,10,10);
DataSens(41:50,5) = linspace(0.1,1,10);
DataSens(51:60,6) = linspace(0.01,0.1,10);
DataSens(61:70,8) = linspace(2,100,10);
DataSens(:,7) = 0.95*DataSens(:,1)./DataSens(:,2).*DataSens(:,8);

R = [];
for ipar = 1:length(I)
    rareload = sprintf('rare_par1000SensAna%d_3_3_2',I(ipar));
    load(rareload);
    R = [R, rare_par];
end

for jpar = 1:70
    L(jpar) = length(find(R == jpar));
end
L = L/length(I);

for iplot = 1:7

figure
if iplot < 7
    plot(DataSens((iplot-1)*10+1:(iplot-1)*10+10,iplot),L((iplot-1)*10+1:(iplot-1)*10+10),'.','Markersize',20,'Color','k')
    hold on
    plot(DataSens((iplot-1)*10+1:(iplot-1)*10+10,iplot),L((iplot-1)*10+1:(iplot-1)*10+10),'-','Linewidth',2,'Color','k')
    hold on
    for j_rarepar = I
        plot([Data1000(j_rarepar,iplot),Data1000(j_rarepar,iplot)],[0,1],'-','Linewidth',2,'Color','r')
        hold on
    end
    xticks([0,max(DataSens((iplot-1)*10+1:(iplot-1)*10+10,iplot))])
else
    plot(DataSens((iplot-1)*10+1:(iplot-1)*10+10,8),L((iplot-1)*10+1:(iplot-1)*10+10),'.','Markersize',20,'Color','k')
    hold on
    plot(DataSens((iplot-1)*10+1:(iplot-1)*10+10,8),L((iplot-1)*10+1:(iplot-1)*10+10),'-','Linewidth',2,'Color','k')
    hold on
    for j_rarepar = I
        plot([Data1000(j_rarepar,8),Data1000(j_rarepar,8)],[0,1],'-','Linewidth',2,'Color','r')
        hold on
    end
    xticks([0,max(DataSens((iplot-1)*10+1:(iplot-1)*10+10,8))])
end
yticks([0,1])
set(0,'DefaultAxesFontName','Arial');
set(gca,'FontSize',12)
box off
set(gca,'linewidth',2)
ylabel('frequency')
xlabel('parameter')
ylim([0,1])
end



