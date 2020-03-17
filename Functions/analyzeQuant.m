%function evaluating several quantitaive traits of simulations with rare
%coordinated high states
%
%INPUT:
%n_species:     number of nodes in the network of given simulation (eg. 2, 3, 5, 8)
%S_outpar:      output of simulation file - number of gene product counts per
%               time unit per gene and state of DNA per gene per time unit (on/off)
%rare_par:      parameter sets giving rise to simulations with rare
%               coordinated high states
%thres:         threshold at which the cell is called to highly express the
%               gene (here: calculated as the 0.8*steady state of high state
%
%OUTPUT:
%
%NumJack:       number of high expressions states per simulation
%AllTimeJack:   total time (units) in which a simulation is in
%               high expression states (\in [1,maxgillespie])
%MaxNumSpecJack:maximum number of simultaniously highly expressed genes per
%               high expression state per simulation 
%TimeJack:      all time points at which the system is in a high expression
%               state
%NumSpecJackReg:total number of genes involved in hugh expression state


%%
function [NumJack, AllTimeJack, NumSpecJackReg, MaxNumSpecJack, TimeJack] = ...
    analyzeQuant(n_species, S_outpar, rare_par, thres)

for iparjack = 1:length(rare_par)
    
    clearvars -except n_species inet isubnet doplot parjack thres iparjack S_outpar rare_par NumJack AllTimeJack NumSpecJackReg MaxNumSpecJack TimeJack
    
    param = rare_par(iparjack);

    f1 = sum(S_outpar{param}(1:n_species,1:end) > thres(param));

    if mod(n_species,2) == 1 %find all time points at which more than half of the species are above the threshold
        L1 = find(f1 >= ceil(n_species/2));
    else
        L1 = find(f1 >= ceil(n_species/2) + 1);
    end
    
    diff_jackpot = diff(L1);
    term_jackpot_ind = find(diff_jackpot > 1); %find jackpot ends (difference greater than 1 between the time points showing jackpot)
    
    L2 = zeros(1,0);
    count_L = 1;
    for idiff = 1:length(term_jackpot_ind)
        if diff_jackpot(term_jackpot_ind(idiff)) < 50 %if time period between the single jackpot regions is below 30, then take it as one jackpot
            if isempty(L2) == 0
                L2(count_L:count_L+diff_jackpot(term_jackpot_ind(idiff))-1) = L1(term_jackpot_ind(idiff))+1:L1(term_jackpot_ind(idiff)+1);
                count_L = length(L2)+1;
            else
                L2(1:count_L+diff_jackpot(term_jackpot_ind(idiff))-1) = L1(term_jackpot_ind(idiff))+1:L1(term_jackpot_ind(idiff)+1);
                count_L = length(L2)+1;
            end
        end
    end
    if isempty(L2) == 0
        L = unique([L1,L2]); %L summarizes all time points at which the system is in a jackpot
    else
        L = L1;
    end
    
    %find jackpot times where times from first time above the threshold to
    %last time above threshold with jackpot in between
    L3 = find(f1 >= 1);
    L3_start = L3([1,find(diff(L3)>50)+1]); %start time of at least one species above threshold
    L3_end = L3([find(diff(L3)>50),length(L3)]); %end time
    
    %find all time point at which at least one species is above threshold
    %and which contains a jackpot (> half the number of nodes is above the
    %threshold) in between L3_start and L3_end
    count_L = 1;
    count_L_start = 1;
    for iL3 = 1:length(L3_start)
        if isempty(intersect(L3_start(iL3):L3_end(iL3),L)) == 0
            L_all(count_L:count_L+length(L3_start(iL3):L3_end(iL3))-1) = L3_start(iL3):L3_end(iL3);
            L_start(count_L_start) = L3_start(iL3);
            L_end(count_L_start) = L3_end(iL3);
            count_L_start = count_L_start + 1;
            count_L = length(L_all) + 1;
        end
    end
    
    time_start_jackpot = L_start;
    time_finish_jackpot = L_end;
    
    %number of jackpots
    NumJack(iparjack) = length(time_start_jackpot);
    
    %overall time in jackpot
    AllTimeJack(iparjack) = length(L_all);
    
    for ijack = 1:length(time_start_jackpot)
        
        %time in each jackpot
        TimeJack{iparjack}(ijack) = time_finish_jackpot(ijack)-time_start_jackpot(ijack);
        
        %number of species involved in a jackpot region
        NumSpecJackReg{iparjack}(ijack) = 0;
        for jspec = 1:n_species
%             if sum(S_outpar{param}(jspec,time_start_jackpot(ijack):time_finish_jackpot(ijack)) > thres(param)) > 1
                if sum(S_outpar{iparjack}(jspec,time_start_jackpot(ijack):time_finish_jackpot(ijack)) > thres(param)) > 1
                NumSpecJackReg{iparjack}(ijack) = NumSpecJackReg{iparjack}(ijack)+1;
            end
        end
        
        %maximum number of species in jackpot at same time
%         MaxNumSpecJack{iparjack}(ijack) = max(sum(S_outpar{param}(1:n_species,time_start_jackpot(ijack):time_finish_jackpot(ijack)) > thres(param)));
        MaxNumSpecJack{iparjack}(ijack) = max(sum(S_outpar{iparjack}(1:n_species,time_start_jackpot(ijack):time_finish_jackpot(ijack)) > thres(param)));
    end
    
end
end