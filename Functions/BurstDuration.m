
load('/Volumes/MELANOMA/Data/Data1000')
thres = Data1000(:,1)./Data1000(:,2).*Data1000(:,8)*0.8;

load('/Volumes/MELANOMAII/Revisions/S_outparPoints3_2_1')
% load('/Volumes/MELANOMAII/Revisions/S_outparPoints1003_2_1')
for n_species = 3
    
    for param = 26
        
        clear base high A B C
        
        %number of species above threshold
        f1 = sum(S_outpar{1}(1:n_species,1:end) > thres(param));
        
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
        
        A = strfind(S_outpar{1}(n_species+2,:),[0,1]); %burst up
        B = strfind(S_outpar{1}(n_species+2,:),[1,0]); %burst down
        
        countbase = 1;
        counthigh = 1;
        for iA = 1:length(A)
            if ismember(A(iA),unique(L)) == 0
                base(countbase) = iA;
                countbase = countbase + 1;
            else
                high(counthigh) = iA;
                counthigh = counthigh + 1;
            end
        end
        
        if length(A) == length(B)
            C = B-A;
        else
            C = B-A(1:end-1);
            base = base(1:end-1);
        end
        Mbase = mean(C(base)-1);
        Mhigh = mean(C(high)-1);
%         hC(countrare) = lillietest(C-1,'Distribution','exponential','Alpha',0.1);
    end
    
end

% figure
% [h,L,MX,MED]=violin([(num_high./lengthnum_high)',(num_base./lengthnum_base)'], 'facecolor', [150,150,150;150,150,150]./255);
% xticks([1,2])
% yticks([0,0.1])
% xticklabels({'high time-region','baseline time-region'})
% set(0,'DefaultAxesFontName','Arial');
% set(gca,'FontSize',12)
% box off
% set(gca,'linewidth',2)
% ylabel('burst frequency')
% ylim([-0.03,0.1])

