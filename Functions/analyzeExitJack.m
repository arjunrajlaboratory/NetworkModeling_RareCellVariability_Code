function sol = analyzeExitJack(n_species, S_outpar, parJ, doplot, inet, thres, isubnet)
if isempty(parJ) == 0
    for iPar = 1:length(parJ)
        
        clearvars -except n_species inet doplot parJ thres iPar S_outpar H P C K pd sol isubnet
        
        param = parJ(iPar);
        
        %number of species above threshold
        f1 = sum(S_outpar{param-(isubnet-1)*length(S_outpar)}(1:n_species,1:end) > thres(param));
        
        %find all time points at which more than half of the species are above
        %the threshold
        if mod(n_species,2) == 1
            L1 = find(f1 >= ceil(n_species/2));
        else
            L1 = find(f1 >= ceil(n_species/2) + 1);
        end
        
        %find jackpot ends (difference greater than 1 between the time points
        %showing jackpot)
        diffJ = diff(L1);
        TermJ_ind = find(diffJ > 1);
        
        L2 = zeros(1,0);
        count_L = 1;
        for iDiff = 1:length(TermJ_ind)
            %if time period between the single jackpot regions is below 50,
            %then take it as one jackpot
            if diffJ(TermJ_ind(iDiff)) < 50
                L2(count_L:count_L+diffJ(TermJ_ind(iDiff))-1) = L1(TermJ_ind(iDiff))+1:L1(TermJ_ind(iDiff)+1);
                count_L = length(L2)+1;
            end
        end
        
        %L summarizes all time points at which the system is in a jackpot region
        if isempty(L2) == 1
            L = unique(L1);
        else
            L = unique([L1,L2]);
        end
        
        %find jackpot times
        L3 = find(f1 >= 1);
        L3_start = L3([1,find(diff(L3)>50)+1]);
        L3_end = L3([find(diff(L3)>50),length(L3)]);
        
        %find all time intervals in which at least one species is above threshold
        %and which contains a jackpot in between start and end of jackpot region
        countL = 1;
        countLStart = 1;
        for iL3 = 1:length(L3_start)
            if isempty(intersect(L3_start(iL3):L3_end(iL3),L)) == 0
                L(countL:countL+length(L3_start(iL3):L3_end(iL3))-1) =...
                    L3_start(iL3):L3_end(iL3);
                LStart(countLStart) = L3_start(iL3);
                LEnd(countLStart) = L3_end(iL3);
                countLStart = countLStart + 1;
                countL = length(L) + 1;
            end
        end
        
        TimeStartJ = LStart;
        TimeFinishJ = LEnd;
        
        DiffTime = NaN(length(TimeFinishJ),n_species-1);
        
        for iTimeFinishJ = 1:length(TimeFinishJ)
            
            %find all species above threshold in jackpot region
            countSpec = 1;
            for jSpec = 1:n_species
                if sum(S_outpar{param-(isubnet-1)*length(S_outpar)}(jSpec,TimeStartJ(iTimeFinishJ):...
                        TimeFinishJ(iTimeFinishJ))>thres(param)) > 0
                    InvSpec{iTimeFinishJ}(countSpec) = jSpec;
                    countSpec = countSpec + 1;
                end
            end
            
            for iSpec = 1:length(InvSpec{iTimeFinishJ})
                LThres{iTimeFinishJ}(iSpec) = TimeStartJ(iTimeFinishJ) ...
                    + find(S_outpar{param-(isubnet-1)*length(S_outpar)}(InvSpec{iTimeFinishJ}(iSpec),...
                    TimeStartJ(iTimeFinishJ):TimeFinishJ(iTimeFinishJ)) > thres(param),...
                    1,'last')-2;
            end
            
            AllDiff{iTimeFinishJ} = diff(sort(LThres{iTimeFinishJ}));
            DiffTime(iTimeFinishJ,1:length(AllDiff{iTimeFinishJ})) = AllDiff{iTimeFinishJ};
        end
        
        sol.Data{iPar} = DiffTime;
        
        %determine for all jackpots
        for iDiffTime = 1:n_species-1
            D = DiffTime(:,iDiffTime);
            D(isnan(D)) = [];
            D1{iDiffTime} = D;
            
            if length(D1{iDiffTime}) > 4
                [h,p,kstat,critval] = lillietest(D1{iDiffTime},'Distr','exp');
                H{iPar}(iDiffTime) = h;
                P{iPar}(iDiffTime) = p;
                K{iPar}(iDiffTime) = kstat;
                C{iPar}(iDiffTime) = critval;
                
                pd{iPar}(iDiffTime) = fitdist(D1{iDiffTime},'Exponential');
            end
        end
        
        if isequal(doplot,'yes')
            f = figure('visible','off');
            for iPlot = 1:n_species-1
                if length(D1{iPlot}) > 4
                    subplot(1,n_species-1,iPlot)
                    if H{iPar}(iPlot) == 0
                        h = histfit(D1{iPlot},[],'exponential');
                        h(1).FaceColor = [0.5 0.5 0.5];
                        h(2).Color = [0 0 0];
                        h(1).BarWidth = 1;
                        xlabel('Time (sec)')
                        ylabel('Count')
                    else
                        histogram(D1{iPlot},'Binwidth', 10, 'Facecolor', [0.8 0.8 0.8])
                    end
                    xlabel('Time (sec)')
                    ylabel('Count')
                end
            end
        end
    end
    if exist('H') == 1
        sol.H = H;
    end
    if exist('P') == 1
        sol.P = P;
    end
    if exist('K') == 1
        sol.K = K;
    end
    if exist('C') == 1
        sol.C = C;
    end
    
    if exist('sol')
%         name_sol = sprintf('/Volumes/MELANOMAII/Example/sol_ExitJack%d_%d_%d',n_species,inet,isubnet);
        name_sol = sprintf('./Example/sol_ExitJack%d_%d_%d',n_species,inet,isubnet);
        save(name_sol,'sol')
    end
    
else
    sol.Data = NaN;
    sol.H = NaN;
    sol.P = NaN;
    sol.K = NaN;
    sol.C = NaN;
%     name_sol = sprintf('/Volumes/MELANOMAII/Example/sol_ExitJack%d_%d_%d',n_species,inet,isubnet);
    name_sol = sprintf('./Example/sol_ExitJack%d_%d_%d',n_species,inet,isubnet);
    save(name_sol,'sol')
end
end