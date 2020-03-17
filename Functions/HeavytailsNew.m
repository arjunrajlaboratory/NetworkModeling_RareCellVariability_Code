
Net = [2,4,10,80];
Subnet = [10,10,100,250];
ncount = 0;

count = 1;
for n_species = [2,3,5,8]
    ncount = ncount +1;
    n_species 
    
    for inet = 1:Net(ncount)
        
%         loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/rare_par1000_%d_%d.mat',n_species,n_species,inet);
        loadrare = sprintf('/Volumes/MELANOMA/Data/RareParameters/%dnodes/Replicates3/rare_par10003_%d_%d.mat',n_species,n_species,inet);
        
        load(loadrare)
        
        for isubnet = 1:Subnet(ncount)
            
            r1 = rare_par(rare_par > isubnet*1000/Subnet(ncount)-1000/Subnet(ncount));
            r2 = rare_par(rare_par <= isubnet*1000/Subnet(ncount));
            rare_par_isubnet = intersect(r1,r2);
            
            if isempty(rare_par_isubnet) == 0
                
%                 loadsol = sprintf('/Volumes/MELANOMAII/Data/CriteriaAnalysis/%dnodes/sol1000%d_%d_%d.mat',n_species,n_species,inet,isubnet);
                loadsol = sprintf('/Volumes/MELANOMAII/Data/CriteriaAnalysis/%dnodes/Replicates3/sol10003%d_%d_%d.mat',n_species,n_species,inet,isubnet);
                load(loadsol)
                
                for i = rare_par_isubnet
                    DataSim = sol{i-(length(sol)*isubnet-length(sol))}.samp(1,:);
                    PrctlSim = prctile(sol{i-(length(sol)*isubnet-length(sol))}.samp(1,:),99);
                    
                    a = histcounts(DataSim,'BinWidth',1);
                    bin = 1;
                    
                    if find(a == max(a)) <= 15
                        
                        anew = mean(DataSim);
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
                        
                    else
                        %0 signals max bins size at value larger 15
                        Com(count,:) = [0,n_species,inet,isubnet,i-(length(sol)*isubnet-length(sol)),i,NaN,NaN];
                        count = count + 1;
                    end
                    
                    
                end
            end
        end
    end
end

A = find(Com(:,7)>Com(:,8));
B = length(find(Com(:,1)==0));
C = sum(isnan(Com(:,7)));

figure
plot(Com(:,8),Com(:,7),'.','MarkerSize',20, 'Color', [211,95,95]./255)
hold on
plot([1,3000],[1,3000],':','Color','k','LineWidth',2)
box off
set(gca,'linewidth',2)
xlabel('99th percentile of fitted exponential distribution')
ylabel({'99th percentile of expression'; 'distributions of simulations with'; 'rare coordinated high states'})
xticks([0,3000])
yticks([0,3000])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)

%% calculate number of simuations per network size

%most stringent
for irep = 1:3
    
    loadCom = sprintf('/Volumes/MELANOMAII/Revisions/Com%d',irep);
    load(loadCom)
    
    count = 1;
    for n_species = [2,3,5,8]
        N(count) = length(find(Com(:,2) == n_species & Com(:,7)>Com(:,8)));
        count = count +1;
    end

    %less stringent
    count = 1;
    N_less1 = find(Com(:,1) == 0);
    for n_species = [2,3,5,8]
        N_less2(count) = length(find(Com(N_less1,2) == n_species));
        count = count +1;
    end
    N_less = N_less2+N;
    
    R(irep,:) = N;
    R_less(irep,:) = N_less;
    
end

figure
boxplot([R(:,1),R(:,2), R(:,3), R(:,4)],'Colors', [0,0,0]./255)
ylim([0,140])
xlabel('network size')
xticklabels({'2','3','5','8'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
ylabel('simulations with rare coordinated high states')
yticks([0,70,140])
box off
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
set(gca,'linewidth',2)

figure
boxplot([R_less(:,1),R_less(:,2), R_less(:,3), R_less(:,4)],'Colors', [0,0,0]./255)
ylim([0,170])
xlabel('network size')
xticklabels({'2','3','5','8'})
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
ylabel('simulations with rare coordinated high states')
yticks([0,85,170])
box off
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
set(gca,'linewidth',2)


%% connectivity for 5 nodes networks

for inet = 1:10
    R(inet) = length(find(Com(:,2) == 5 & Com(:,3) == inet & Com(:,7)>Com(:,8)));
end

load('/Volumes/MELANOMA/Data/M_iso5')
for jnet = 1:length(M_iso)
    C(jnet) = sum(M_iso{jnet}(:,1));
end

figure
plot(C,R,'.','Markersize',20,'Color', [211,95,95]./255)
ylabel('simulations with rare coordinated high states')
xlabel('connectivity')
xticks([1,2,3,4,5])
yticks([0,20,30])
xlim([0.5,5.5])
ylim([0,30])
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca,'linewidth',2)

%less stringent

N_less1 = find(Com(:,1) == 0);
for inet = 1:10
    R(inet) = length(find(Com(:,2) == 5 & Com(:,3) == inet & Com(:,7)>Com(:,8)));
    R(inet) = R(inet)+length(find(Com(N_less1,2) == 5 & Com(N_less1,3) == inet));
end

load('/Volumes/MELANOMA/Data/M_iso5')
for jnet = 1:length(M_iso)
    C(jnet) = sum(M_iso{jnet}(:,1));
end

figure
plot(C,R,'.','Markersize',20,'Color', [211,95,95]./255)
ylabel('simulations with rare coordinated high states')
xlabel('connectivity')
xticks([1,2,3,4,5])
yticks([0,20,40])
xlim([0.5,5.5])
ylim([0,40])
box off
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
set(gca,'linewidth',2)

%% jackpot parameter sets

A = find(Com(:,7)>Com(:,8));
for inum = 1:1000
    ParFreq(inum) = length(find(Com(A,6) == inum));
end

figure
histogram(ParFreq(ParFreq<15)/96*100,'FaceColor',[0.2, 0.2, 0.2],'Binwidth',1,...
    'Edgecolor','none')
hold on
histogram(ParFreq(ParFreq>=15)/96*100,'FaceColor',[255, 127, 42]./255','Edgecolor','none','Binwidth',1)
alpha(1)
hold on 
plot([15,15],[0,970],'--','Color','k','Linewidth',3)
xlim([0,50])
ylim([0,970])
xlabel('% of simulations with rare coordinated high states per parameter set')
ylabel('Counts')
set(gca,'linewidth',2)
box off
xticks([0,25,50])
yticks([0, 10, 20, 30, 960, 970])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
breakyaxis([30, 960]);

%less stringent
N_less1 = find(Com(:,1) == 0);
A = find(Com(:,7)>Com(:,8));
for inum = 1:1000
    ParFreq(inum) = length(find(Com(A,6) == inum));
    ParFreq(inum) = ParFreq(inum) + length(find(Com(N_less1,1) == 0 & Com(N_less1,6) == inum));
end

figure
histogram(ParFreq(ParFreq<15)/96*100,'FaceColor',[0.2, 0.2, 0.2],'Binwidth',1,...
    'Edgecolor','none')
hold on
histogram(ParFreq(ParFreq>=15)/96*100,'FaceColor',[255, 127, 42]./255','Edgecolor','none','Binwidth',1)
alpha(1)
hold on 
plot([15,15],[0,920],'--','Color','k','Linewidth',3)
xlim([0,50])
ylim([0,940])
xlabel('% of simulations with rare coordinated high states per parameter set')
ylabel('Counts')
set(gca,'linewidth',2)
box off
xticks([0,25,50])
yticks([0, 10, 20, 30, 930, 940])
set(0,'DefaultAxesFontName','Arial'); 
set(gca,'FontSize',12)
breakyaxis([30, 930]);
