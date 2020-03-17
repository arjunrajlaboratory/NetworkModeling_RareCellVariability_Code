function analyzeBurstsTerm(Nspecies, S_outpar, rare_par, inet, thres, isubnet)

RareCount = rare_par;
if isempty(RareCount) == 0   
for iRareCount = 1:length(RareCount)
    
    clearvars -except Nspecies S_outpar RareCount inet thres iRareCount...
        BurstOffJExit BurstOffJNExit BurstOnJExit BurstOnJNExit...
        BurstOffJExit_Data BurstOffJNExit_Data BurstOnJExit_Data...
        BurstOnJNExit_Data isubnet
    
    %parameter set to be tested
    param = RareCount(iRareCount);
    
    %determine the number of species above thres
    f1 = sum(S_outpar{param-(isubnet-1)*length(S_outpar)}(1:Nspecies,1:end) > thres(param));
    
    %find jackpot time points
    %jackpot ==(more than half of the species must be above the threshold)
    if mod(Nspecies,2) == 1
        L1 = find(f1 >= ceil(Nspecies/2));
    else
        L1 = find(f1 >= ceil(Nspecies/2) + 1);
    end
    
    %if the gap between two jackpots is < than 50, consider them to
    %belong together (jackpot region)
    L2 = zeros(0,1);
    diffJ = diff(L1);
    TermJ_ind = find(diffJ > 1);
    countL = 1;
    for iDiff = 1:length(TermJ_ind)
        if diffJ(TermJ_ind(iDiff)) < 50
            L2(1,countL:countL+diffJ(TermJ_ind(iDiff))-1)...
                = L1(TermJ_ind(iDiff))+1:L1(TermJ_ind(iDiff)+1);
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
    L3_start = L3([1,find(diff(L3)>50)+1]);
    L3_end = L3([find(diff(L3)>50),length(L3)]);
    
    %find all time intervals in which at least one species is above threshold
    %and which contains a jackpot in between start and end of jackpot region
    countL = 1;
    countL_start = 1;
    for iL3 = 1:length(L3_start)
        if isempty(intersect(L3_start(iL3):L3_end(iL3),L)) == 0
            LAll(countL:countL+length(L3_start(iL3):L3_end(iL3))-1) =...
                L3_start(iL3):L3_end(iL3);
            LStart(countL_start) = L3_start(iL3);
            LEnd(countL_start) = L3_end(iL3);
            countL_start = countL_start + 1;
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
    
    %determine exiting region for each jackpot region
    %find first species exiting (and staying below) gene high state within
    %last quarter of jackpot region
    for iJ = 1:length(TimeStartJ)
        
        clear LastHigh
        
        %find all species above threshold in last quarter of jackpot region
        SpecExit{iJ} = find(sum(S_outpar{param-(isubnet-1)*length(S_outpar)}(1:Nspecies,max(1,...
            TimeFinishJ(iJ)-ceil((TimeFinishJ(iJ)-TimeStartJ(iJ)+1)/4)):...
            TimeFinishJ(iJ)) > thres(param),2) > 0);
        
        %find fist species leaving high-expression state and not returning
        for iSpec = 1:length(SpecExit{iJ})
            LastHigh(iSpec) = TimeFinishJ(iJ)-ceil((TimeFinishJ(iJ)-...
                TimeStartJ(iJ)+1)/4) + find(S_outpar{param-(isubnet-1)*length(S_outpar)}(SpecExit{iJ}(iSpec),...
                max(1,TimeFinishJ(iJ)-ceil((TimeFinishJ(iJ)-TimeStartJ(iJ)+1)/4)):...
                TimeFinishJ(iJ)) > thres(param),1,'last')-2;
        end
        BeginExit(iJ) = min(LastHigh);
    end
    
    for iSpec = 1:Nspecies %consider all bursts for each species seperate
        
        clearvars -except Nspecies S_outpar RareCount inet thres iRareCount...
            BurstOffJExit BurstOffJNExit BurstOnJExit BurstOnJNExit...
            BurstOffJExit_Data BurstOffJNExit_Data BurstOnJExit_Data...
            BurstOnJNExit_Data LAll param LAll TimeStartJ TimeFinishJ iSpec...
            BeginExit isubnet
        
        countOff = 0;
        countOn = 0;
        for jJ = 1:length(TimeStartJ)
            
            %find all burst off states in jackpot region
            BurstOffJ = find(S_outpar{param-(isubnet-1)*length(S_outpar)}(Nspecies-1+2*iSpec,...
                TimeStartJ(jJ):TimeFinishJ(jJ)) == 1);
            for iOff = 1:length(BurstOffJ)
                if iOff == 1
                    countOff = countOff + 1;
                    countOffNew = countOff;
                    BurstOffTimesJ{countOff} = TimeStartJ(jJ)-1+BurstOffJ(iOff);
                elseif BurstOffJ(iOff)-1 == BurstOffJ(iOff-1)
                    BurstOffTimesJ{countOff} = [BurstOffTimesJ{countOff},...
                        TimeStartJ(jJ)-1+BurstOffJ(iOff)];
                else
                    countOff = countOff +1;
                    BurstOffTimesJ{countOff} = TimeStartJ(jJ)-1+BurstOffJ(iOff);
                end
            end
            
            if exist('countOffNew') == 1
                for iCountOff = countOffNew:countOff
                    BurstOffTimeLengthsJ(iCountOff) = length(BurstOffTimesJ{iCountOff});
                    if isempty(intersect(BurstOffTimesJ{iCountOff},...
                            BeginExit(jJ):TimeFinishJ(jJ))) == 0
                        if isempty(intersect(BurstOffTimesJ{iCountOff}, TimeStartJ)) == 1
                            BurstOffExit(iCountOff) = 1;
%                         elseif isempty(intersect(BurstOffTimesJ{iCountOff}, TimeFinishJ)) == 1
%                             BurstOffExit(iCountOff) = 1;
                        else
                            BurstOffExit(iCountOff) = 2;
                        end
                    else
                        if isempty(intersect(BurstOffTimesJ{iCountOff}, TimeStartJ)) == 1
                            BurstOffExit(iCountOff) = 0;
%                         elseif isempty(intersect(BurstOffTimesJ{iCountOff}, TimeFinishJ)) == 1
%                             BurstOffExit(iCountOff) = 0;
                        else
                            BurstOffExit(iCountOff) = 2;
                        end
                    end
                end
            end
            
            %find all burst on states in jackpot region
            BurstOnJ = find(S_outpar{param-(isubnet-1)*length(S_outpar)}(Nspecies+2*iSpec,...
                TimeStartJ(jJ):TimeFinishJ(jJ)) == 1);
            for iOn = 1:length(BurstOnJ)
                if iOn == 1
                    countOn = countOn + 1;
                    countOnNew = countOn;
                    BurstOnTimesJ{countOn} = TimeStartJ(jJ)-1+BurstOnJ(iOn);
                elseif BurstOnJ(iOn)-1 == BurstOnJ(iOn-1)
                    BurstOnTimesJ{countOn} = [BurstOnTimesJ{countOn},...
                        TimeStartJ(jJ)-1+BurstOnJ(iOn)];
                else
                    countOn = countOn +1;
                    BurstOnTimesJ{countOn} = TimeStartJ(jJ)-1+BurstOnJ(iOn);
                end
            end
            
            if exist('countOnNew') == 1
                for iCountOn = countOnNew:countOn
                    BurstOnTimeLengthsJ(iCountOn) = length(BurstOnTimesJ{iCountOn});
                    if isempty(intersect(BurstOnTimesJ{iCountOn},...
                            BeginExit(jJ):TimeFinishJ(jJ))) == 0
                        if isempty(intersect(BurstOnTimesJ{iCountOn}, TimeStartJ)) == 1
                            BurstOnExit(iCountOn) = 1;
%                         elseif isempty(intersect(BurstOnTimesJ{iCountOn}, TimeFinishJ)) == 1
%                             BurstOnExit(iCountOn) = 1;
                        else
                            BurstOnExit(iCountOn) = 2;
                        end
                    else
                        if isempty(intersect(BurstOnTimesJ{iCountOn}, TimeStartJ)) == 1
                            BurstOnExit(iCountOn) = 0;
%                         elseif isempty(intersect(BurstOnTimesJ{iCountOn}, TimeFinishJ)) == 1
%                             BurstOnExit(iCountOn) = 0;
                        else
                            BurstOnExit(iCountOn) = 2;
                        end
                    end
                end
            end
        end
        
        if exist('BurstOffExit') == 1
            ExitBurstOffs = find(BurstOffExit == 1);
            BurstOffJExit_Data{iRareCount,iSpec} = BurstOffTimeLengthsJ(ExitBurstOffs);
            BurstOffJExit(iRareCount,iSpec) = median(BurstOffTimeLengthsJ(ExitBurstOffs));
            
            NExitBurstOffs = find(BurstOffExit == 0);
            BurstOffJNExit_Data{iRareCount,iSpec} = BurstOffTimeLengthsJ(NExitBurstOffs);
            BurstOffJNExit(iRareCount,iSpec) = median(BurstOffTimeLengthsJ(NExitBurstOffs));
        else
            BurstOffJExit_Data{iRareCount,iSpec} = [];
            BurstOffJExit(iRareCount,iSpec) = NaN;
            
            BurstOffJNExit_Data{iRareCount,iSpec} = [];
            BurstOffJNExit(iRareCount,iSpec) = NaN;
        end
        
        if exist('BurstOnExit') == 1
            ExitBurstOns = find(BurstOnExit == 1);
            BurstOnJExit_Data{iRareCount,iSpec} = BurstOnTimeLengthsJ(ExitBurstOns);
            BurstOnJExit(iRareCount,iSpec) = median(BurstOnTimeLengthsJ(ExitBurstOns));
            
            NExitBurstOns = find(BurstOnExit == 0);
            BurstOnJNExit_Data{iRareCount,iSpec} = BurstOnTimeLengthsJ(NExitBurstOns);
            BurstOnJNExit(iRareCount,iSpec) = median(BurstOnTimeLengthsJ(NExitBurstOns));
        else
            BurstOnJExit_Data{iRareCount,iSpec} = [];
            BurstOnJExit(iRareCount,iSpec) = NaN;
            
            BurstOnJNExit_Data{iRareCount,iSpec} = [];
            BurstOnJNExit(iRareCount,iSpec) = NaN;
        end
    end
end

%save all
M.BurstOffJExit_Data = BurstOffJExit_Data;
M.BurstOffJNExit_Data = BurstOffJNExit_Data;
M.BurstOnJExit_Data = BurstOnJExit_Data;
M.BurstOnJNExit_Data = BurstOnJNExit_Data;


for iStat = 1:size(BurstOffJExit,1)
    %Mann-Whitney U-test/Wilcoxon rank tum test
    if length([BurstOffJExit_Data{iStat,:}]) > 10
        [p,h] = ranksum([BurstOffJExit_Data{iStat,:}],[BurstOffJNExit_Data{iStat,:}]);
        StatSigOff(iStat) = double(h);
        POff(iStat) = p;
    else
        StatSigOff(iStat) = NaN;
        POff(iStat) = NaN;
    end
    %ks2statoff(istat) = ks2stat;
    MedianBurstOffJExit(iStat) = median([BurstOffJExit_Data{iStat,:}]);
    MedianBurstOffJNExit(iStat) = median([BurstOffJNExit_Data{iStat,:}]);
end

for iStat = 1:size(BurstOnJExit,1)
    if length([BurstOnJExit_Data{iStat,:}]) > 1
        %Kolmogorov-Smirnoff test
        [h,p,ks2stat] = kstest2([BurstOnJExit_Data{iStat,:}],...
            [BurstOnJNExit_Data{iStat,:}]);
        StatSigOn(iStat) = double(h);
        POn(iStat) = p;
        ks2StatOn(iStat) = ks2stat;
    else
        StatSigOn(iStat) = NaN;
        POn(iStat) = NaN;
        ks2StatOn(iStat) = NaN;
    end
end

M.StatSigOn = StatSigOn';
M.StatSigOff = StatSigOff';
M.POn = POn';
M.POff = POff';
M.ks2StatOn = ks2StatOn';

M.MedianBurstOffJExit = MedianBurstOffJExit';
M.MedianBurstOffJNExit = MedianBurstOffJNExit';
M.Par = RareCount;

% save_MTerm = sprintf('/Volumes/MELANOMAII/Example/sol_BurstsTerm%d_%d_%d',...
%     Nspecies,inet, isubnet);
save_MTerm = sprintf('./Example/sol_BurstsTerm%d_%d_%d',...
    Nspecies,inet, isubnet);
save(save_MTerm,'M')

else
    M.BurstOffJExit_Data = NaN;
    M.BurstOffJNExit_Data = NaN;
    M.BurstOnJExit_Data = NaN;
    M.BurstOnJNExit_Data = NaN;
    M.StatSigOn = NaN;
    M.StatSigOff = NaN;
    M.POn = NaN;
    M.POff = NaN;
    M.ks2StatOn = NaN;

    M.MedianBurstOffJExit = NaN;
    M.MedianBurstOffJNExit = NaN;
    M.Par = NaN;
%     save_MTerm = sprintf('/Volumes/MELANOMAII/Example/sol_BurstsTerm%d_%d_%d',...
%     Nspecies,inet, isubnet);
    save_MTerm = sprintf('./Example/sol_BurstsTerm%d_%d_%d',...
    Nspecies,inet, isubnet);
    save(save_MTerm,'M')

end
end
