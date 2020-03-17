function analyzeBurstsInit(Nspecies, S_outpar, rare_par, inet, thres, isubnet)

RareCount = rare_par;
if isempty(RareCount) == 0
    for iRareCount = 1:length(RareCount)
        
        clearvars -except  Nspecies S_outpar RareCount inet thres iRareCount...
            BurstOnJInit BurstOnJNInit BurstOnJInit_Data BurstOnJNInit_Data isubnet
        
        %parameter set to be tested
        param = RareCount(iRareCount);
        
        %determine the number of species above thres
        f1 = sum(S_outpar{param-(isubnet-1)*length(S_outpar)}(1:Nspecies,1:end) > thres(param));
        
        %find jackpot times
        %jackpot == more than half of the species must be above the threshold
        if mod(Nspecies,2) == 1
            L1 = find(f1 >= ceil(Nspecies/2));
        else
            L1 = find(f1 >= ceil(Nspecies/2) + 1);
        end
        
        %if the gap between two jackpots is < than 50, consider them to
        %belong together (jackpot region)
        L2 = zeros(0,1);
        DiffJ = diff(L1);
        TermJ_ind = find(DiffJ > 1);
        countL = 1;
        for idiff = 1:length(TermJ_ind)
            if DiffJ(TermJ_ind(idiff)) < 50
                L2(1,countL:countL+DiffJ(TermJ_ind(idiff))-1)...
                    = L1(TermJ_ind(idiff))+1:L1(TermJ_ind(idiff)+1);
                countL = length(L2)+1;
            end
        end
        
        %L summarizes all time points at which the system is in a jackpot region
        if isempty(L2) == 1
            L = unique(L1);
        else
            L = unique([L1,L2]);
        end
        
        %find all times at which we exceed threshold with at least one species:
        %again 50 sec gaps in between are considered to be one jackpot region
        L3 = find(f1 >= 1);
        L3Start = L3([1,find(diff(L3)>50)+1]);
        L3End = L3([find(diff(L3)>50),length(L3)]);
        
        %find all time intervals in which at least one species is above threshold
        %and which contains a jackpot in between start and end of jackpot region
        countL = 1;
        countLStart = 1;
        for iL3 = 1:length(L3Start)
            if isempty(intersect(L3Start(iL3):L3End(iL3),L)) == 0
                LAll(countL:countL+length(L3Start(iL3):L3End(iL3))-1) =...
                    L3Start(iL3):L3End(iL3);
                LStart(countLStart) = L3Start(iL3);
                LEnd(countLStart) = L3End(iL3);
                countLStart = countLStart + 1;
                countL = length(LAll) + 1;
            end
        end
        
        if LStart(1) == 1
            TimeStartJ = LStart(2:end);
            TimeFinishJ = LEnd(2:end);
        else
            TimeStartJ = LStart;
            TimeFinishJ = LEnd;
        end
        
        %all jackpot region time points
        for iTimesJ = 1:length(TimeStartJ)
            if iTimesJ == 1
                TimesJ(1:TimeFinishJ(iTimesJ)-TimeStartJ(iTimesJ)+1)  = ...
                    TimeStartJ(iTimesJ):TimeFinishJ(iTimesJ);
            else
                TimesJ(length(TimesJ)+1:length(TimesJ)+...
                    TimeFinishJ(iTimesJ)-...
                    TimeStartJ(iTimesJ)+1) = ...
                    TimeStartJ(iTimesJ):TimeFinishJ(iTimesJ);
            end
        end
        
        %find the species initiating the jackpot region
        for iLstart = 1:length(TimeStartJ)
            if TimeStartJ(iLstart) ~= 1
                SpecJ_ind(iLstart) =...
                    find(S_outpar{param-(isubnet-1)*length(S_outpar)}(1:Nspecies,TimeStartJ(iLstart))...
                    > thres(param),1);
            end
        end
        
        for iSpec = 1:Nspecies %consider all bursts for each species seperate
            
            if exist('TimesJ')
            
                clearvars -except Nspecies S_outpar RareCount inet thres iRareCount...
                    BurstOnJInit BurstOnJNInit BurstOnJInit_Data BurstOnJNInit_Data...
                    param TimeStartJ TimeFinishJ iSpec TimesJ SpecJ_ind isubnet

                CountOn = 0;
                for jTime = 1:length(S_outpar{param-(isubnet-1)*length(S_outpar)})

                    %find all burst-on states not in jackpot region
                    if S_outpar{param-(isubnet-1)*length(S_outpar)}(Nspecies+2*iSpec,jTime) == 1
                        %begin burst
                        if S_outpar{param-(isubnet-1)*length(S_outpar)}(Nspecies+2*iSpec,jTime-1) == 0
                            %do not begin in jackpot region
                            if isempty(intersect(jTime, TimesJ)) == 1
                                CountOn = CountOn + 1;
                                BurstOnTimesJ{CountOn}(1) = jTime;
                            end
                        end
                        %continue burst
                        if exist('BurstOnTimesJ') == 1
                            if S_outpar{param-(isubnet-1)*length(S_outpar)}(Nspecies+2*iSpec,jTime-1) == 1
                                BurstOnTimesJ{CountOn}(length(BurstOnTimesJ{CountOn})+1) = jTime;
                            end
                        end
                    end
                end

                %split bursts into 'jackpot initiating' and 'not jackpot
                %initiating' bursts
                for iCount = 1:length(BurstOnTimesJ)
                    StartBurstOn(iCount) = BurstOnTimesJ{iCount}(1);
                    BurstOnTimeLengths(iCount) = length(BurstOnTimesJ{iCount});
                end

                %find indices of the jackpots initiated by this species (ispec)
                InitJSpec_ind = find(SpecJ_ind == iSpec);
                %find the time points at which the species exceeds the threshold
                %the first time before full jackpot
                InitJTimeSpec = TimeStartJ(InitJSpec_ind);

                %find the corresponding burst
                for iJ = 1:length(InitJTimeSpec)
                    diffTime = InitJTimeSpec(iJ)-StartBurstOn;
                    indInitBurst(iJ) = find(InitJTimeSpec(iJ) - min(diffTime(diffTime>0))...
                        == StartBurstOn);
                end

                if exist('indInitBurst')
                    indNInitBurst = setdiff(1:length(BurstOnTimesJ),indInitBurst);
                    BurstOnJInit_Data{iRareCount,iSpec} = BurstOnTimeLengths(indInitBurst);
                    BurstOnJInit(iRareCount,iSpec) = median(BurstOnTimeLengths(indInitBurst));
                    BurstOnJNInit_Data{iRareCount,iSpec} = BurstOnTimeLengths(indNInitBurst);
                    BurstOnJNinit(iRareCount,iSpec) = median(BurstOnTimeLengths(indNInitBurst));

                else
                    indNInitBurst = 1:length(BurstOnTimesJ);
                    BurstOnJInit_Data{iRareCount,iSpec} = [];
                    BurstOnJInit(iRareCount,iSpec) = NaN;
                    BurstOnJNInit_Data{iRareCount,iSpec} = BurstOnTimeLengths(indNInitBurst);
                    BurstOnJNinit(iRareCount,iSpec) = median(BurstOnTimeLengths(indNInitBurst));
                end
            else
                indNInitBurst = NaN;
                BurstOnJInit_Data{iRareCount,iSpec} = [];
                BurstOnJInit(iRareCount,iSpec) = NaN;
                BurstOnJNInit_Data{iRareCount,iSpec} = NaN;
                BurstOnJNinit(iRareCount,iSpec) = NaN;
            end
        end
    end
    
    %save all
    M.BurstOnJInit = BurstOnJInit;
    M.BurstOnJInit_Data = BurstOnJInit_Data;
    M.BurstOnJNInit =  BurstOnJNinit;
    M.BurstOnJNInit_Data = BurstOnJNInit_Data;
    
    if size(BurstOnJInit,1) > 0
        for iStat = 1:size(BurstOnJInit,1)
            if isempty(BurstOnJInit_Data{iStat,1}) == 0
                [h,p,ks2stat] = kstest2([BurstOnJInit_Data{iStat,:}],...
                    [BurstOnJNInit_Data{iStat,:}]);
                StatSigOn(iStat) = h;
                POn(iStat) = p;
                ks2StatOn(iStat) = ks2stat;
            end
        end

        M.StatSig = StatSigOn';
        M.P = POn';
        M.ks2Stat = ks2StatOn';

        M.Par = RareCount;

%         save_MInit = sprintf('/Volumes/MELANOMAII/Example/sol_BurstsOnInit%d_%d_%d',Nspecies,inet,isubnet);
        save_MInit = sprintf('./Example/sol_BurstsOnInit%d_%d_%d',Nspecies,inet,isubnet);
        save(save_MInit,'M')
    else
        M.BurstOnJInit = NaN;
        M.BurstOnJInit_Data = NaN;
        M.BurstOnJNInit =  NaN;
        M.BurstOnJNInit_Data = NaN;
%         save_MInit = sprintf('/Volumes/MELANOMAII/Example/sol_BurstsOnInit%d_%d_%d',Nspecies,inet,isubnet);
        save_MInit = sprintf('./Example/sol_BurstsOnInit%d_%d_%d',Nspecies,inet,isubnet);
        save(save_MInit,'M')
    end
else
    M.BurstOnJInit = NaN;
    M.BurstOnJInit_Data = NaN;
    M.BurstOnJNInit =  NaN;
    M.BurstOnJNInit_Data = NaN;
%     save_MInit = sprintf('/Volumes/MELANOMAII/Example/sol_BurstsOnInit%d_%d_%d',Nspecies,inet,isubnet);
    save_MInit = sprintf('./Example/sol_BurstsOnInit%d_%d_%d',Nspecies,inet,isubnet);
    save(save_MInit,'M')
end
end
