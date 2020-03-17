%example of supplementary analysis - for 50 parameter sets

% if ~exist('/Volumes/MELANOMAII/Example', 'dir')
%     mkdir /Volumes/MELANOMAII/Example
% end
if ~exist('./Example', 'dir')
    mkdir ./Example
end
addpath(genpath(pwd))

addpath(genpath('/Volumes/MELANOMAII'))

%% Generate Data - Old Model (SupplementaryFigure1)

%%%%%%%%%FILL OUT%%%%%%%%%%%
nruns = 50;                             %corresponds to number of parameter sets to be run
n_species = 2;                          %number of nodes
maxgillespie = 1000000;                 %number of Gillespie simulated time units
gen = 'no';                             %yes = generate new parameters (or 'no' = load parameters)
init.spec = repmat(20,1,n_species);     %initial values for all species (at t=0)
init.Bon = zeros(1,n_species);          %initial values for the burst (where 0 = 'off' and 1 = 'on')
data = 'Data50';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

generateDataModelOld(nruns,n_species,maxgillespie,gen,init,data)

%% Heavytails (SupplementaryFigure2)

clearvars
clc

%%%%%%%%FILL OUT%%%%%%%%%%%
n_species = 2;
Net = 2;
Subnet = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 1;

for inet = 1:Net
    %     loadrare = sprintf('/Volumes/MELANOMAII/Example/rare_par%d_%d.mat',...
    %         n_species,inet);
    loadrare = sprintf('./Example/rare_par%d_%d.mat',...
        n_species,inet);
    load(loadrare)
    
    for isubnet = 1:Subnet
        
        r1 = rare_par(rare_par > isubnet*1000/Subnet-1000/Subnet);
        r2 = rare_par(rare_par <= isubnet*1000/Subnet);
        rare_par_isubnet = intersect(r1,r2);
        
        if isempty(rare_par_isubnet) == 0
            
            %             loadsol = sprintf('Volumes/MELANOMAII/Example/sol%d_%d_%d.mat',...
            %                 n_species,inet,isubnet);
            loadsol = sprintf('./Example/sol%d_%d_%d.mat',...
                n_species,inet,isubnet);
            load(loadsol)
            
            for i = rare_par_isubnet
                DataSim = sol{i-(length(sol)*isubnet-length(sol))}.samp(1,:);
                
                PrctlSim = prctile(sol{i-(length(sol)*isubnet-length(sol))}.samp(1,:),99);
                
                a = histcounts(DataSim,'BinWidth',1);
                bin = 1;
                
                if find(a == max(a)) <= 15
                    
                    anew = max(a);
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
                    Com(count,:) = [0,n_species,inet,isubnet,i-(length(sol)*isubnet-length(sol)),i,NaN,NaN];  %NaN for max bin > 15
                    count = count+1;
                    
                    anew = max(a);
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
                    
                end
                
            end
        end
    end
end

% save('/Volumes/MELANOMAII/Example/Com','Com')
save('./Example/Com','Com')

%% Quantitative Analysis (SupplementaryFigure2)

clearvars;
clc;

%%%%%%%%%FILL OUT%%%%%%%%%%%
n_species = 2;
Net = 2;
Subnet = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for inet = 1:Net
    count = 1;
    rare_par = zeros(1,0);
    S = zeros(1,0);
    
    for isubnet = 1:Subnet
        
        %         loadS = sprintf('/Volumes/MELANOMAII/Example/S_outpar%d_%d_%d',...
        %             n_species,inet,isubnet);
        loadS = sprintf('./Example/S_outpar%d_%d_%d',...
            n_species,inet,isubnet);
        %         loadrare = sprintf('/Volumes/MELANOMAII/Example/rare_par%d_%d',...
        %             n_species,inet);
        loadrare = sprintf('./Example/rare_par%d_%d',...
            n_species,inet);
        
        load(loadS);
        load(loadrare);
        
    end
    
    if isempty(rare_par) == 0
        
        %         load('/Volumes/MELANOMAII/Example/Data50')
        load('./Example/Data50')
        thres = Data50(:,1)./Data50(:,2).*Data50(:,8)*0.8;
        
        [NumJack_sol, AllTimeJack_sol, NumSpecJackReg_sol, MaxNumSpecJack_sol, TimeJack_sol] = ...
            analyzeQuant(n_species, S_outpar, rare_par, thres);
        
        solQuant.NumJack = NumJack_sol;
        solQuant.AllTimeJack = AllTimeJack_sol;
        solQuant.NumSpecJackReg = NumSpecJackReg_sol;
        solQuant.MaxNumSpecJack = MaxNumSpecJack_sol;
        solQuant.TimeJack = TimeJack_sol;
        
    else
        solQuant.NumJack = zeros(1,0);
        solQuant.AllTimeJack = zeros(1,0);
        solQuant.NumSpecJackReg = zeros(1,0);
        solQuant.MaxNumSpecJack = zeros(1,0);
        solQuant.TimeJack = zeros(1,0);
    end
    
    %     save_sol = sprintf('/Volumes/MELANOMAII/Example/solQuant1000%d_%d',...
    %         n_species,inet);
    save_sol = sprintf('./Example/solQuant1000%d_%d',...
        n_species,inet);
    save(save_sol,'solQuant');
end

%% Generate Data - Asymmetric Analysis - Architecture (SupplementaryFigure4)

clearvars
clc

%asymmetric architecture

%%%%%%%%%FILL OUT%%%%%%%%%%%
nruns = 50;                            %corresponds to number of parameter sets to be run
n_species = 2;                          %number of nodes
maxgillespie = 1000000;                 %number of Gillespie simulated seconds
gen = 'no';                             %yes = generate new parameters (or 'no' = load parameters)
init.spec = repmat(20,1,n_species);     %initial values for all species (at t=0)
init.Bon = zeros(1,n_species);          %initial values for the burst (where 0 = 'off' and 1 = 'on')
data = 'Data50';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%generate random network
%how many edges in network
NumEdges = randi([1,n_species^2],1,1);
Edges = randsample(n_species^2,NumEdges);
Mat = zeros(n_species);
Mat(Edges) = 1;

%check if network is connected
G = digraph(Mat);
bins = conncomp(G,'Type','weak');

% save('/Volumes/MELANOMAII/Example/AsymArchMat','Mat')
save('./Example/AsymArchMat','Mat')

generateData(nruns,n_species,maxgillespie,gen,init,data,'asym')

%% Population level - Asymmetric Architecture (SupplementaryFigure4)

clearvars
clc

%%%%%%%%FILL OUT%%%%%%%%%%%
n_species = 2;
Subnet = 1;
inet = 1;
doplot = 'no';
tcell = 1000;                          %time per 'cell'
%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 1;

% load('/Volumes/MELANOMAII/Example/Data50')
load('./Example/Data50')
thres = Data50(:,1)./Data50(:,2).*Data50(:,8)*0.8;

for isubnet = 1:Subnet
    
    clear sol
    
    %         loadS = sprintf('/Volumes/MELANOMAII/Example/S_outparAsymArch%d_%d_%d',...
    %             n_species,inet,isubnet);
    loadS = sprintf('./Example/S_outparAsymArch%d_%d_%d',...
        n_species,inet,isubnet);
    load(loadS);
    
    for i = 1:length(S_outpar)
        
        param = (isubnet-1)*length(S_outpar)+i;
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
    
    for j = 1:length(S_outpar)
        
        if  sol{j}.maxjackpot == 1
            if sol{j}.desc == 1
                if sol{j}.rightskewed == 1
                    if sol{j}.unimodal == 1
                        rare_par(count) = (isubnet-1)*10+j;
                        count = count + 1;
                    end
                end
            end
        end
    end
    
    %         save_sol = sprintf('/Volumes/MELANOMAII/Example/solAsymArch%d_%d_%d',n_species,inet,isubnet);
    save_sol = sprintf('./Example/solAsymArch%d_%d_%d',n_species,inet,isubnet);
    save(save_sol,'sol');
    
end

if exist('rare_par')
    %         saverarepar = sprintf('/Volumes/MELANOMAII/Example/rare_parAsymArch%d_%d',n_species,inet);
    saverarepar = sprintf('./Example/rare_parAsymArch%d_%d',n_species,inet);
    save(saverarepar,'rare_par')
    rare_par
else
    rare_par = [];
    %         saverarepar = sprintf('/Volumes/MELANOMAII/Example/rare_parAsymArch%d_%d',n_species,inet);
    saverarepar = sprintf('./Example/rare_parAsymArch%d_%d',n_species,inet);
    save(saverarepar,'rare_par')
    'No rare_par exists'
end

%% Generate Data - Asymmetric Analysis - Parameter Set (SupplementaryFigure4)

clearvars
clc

%%%%%%%%%FILL OUT%%%%%%%%%%%
n_runs = 50;
n_species = 2;                     %number of nodes
maxgillespie = 1000000;                 %number of Gillespie simulated seconds
gen = 'no';                             %yes = generate new parameters (or 'no' = load parameters)
init.spec = repmat(20,1,n_species);     %initial values for all species (at t=0)
init.Bon = zeros(1,n_species);          %initial values for the burst (where 0 = 'off' and 1 = 'on')
data = 'Data50Asym';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sample 10.000 parameter sets for a f node network
Data50Asym = LHSsamplingAsymPar(n_runs,n_species);

% save('/Volumes/MELANOMAII/Example/Data50Asym','Data50Asym')
save('./Example/Data50Asym','Data50Asym')

generateDataAsymPar(n_runs,n_species,maxgillespie,gen,init,data);

%% Population level - Asymmetric Parameter Set (SupplementaryFigure4)

clearvars
clc

%%%%%%%%FILL OUT%%%%%%%%%%%
n_species = 2;
Subnet = 1;
inet = 1;
doplot = 'no';
tcell = 1000;                          %time per 'cell'
%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 1;

% load('/Volumes/MELANOMAII/Example/Data50Asym')
load('./Example/Data50Asym')
thres = Data50Asym(:,1:n_species)./Data50Asym(:,n_species+1:2*n_species).*Data50Asym(:,7*n_species+1:8*n_species)*0.8;

for isubnet = 1:Subnet
    
    clear sol
    
    %         loadS = sprintf('/Volumes/MELANOMAII/Example/S_outparAsymPar%d_%d_%d',...
    %             n_species,inet,isubnet);
    loadS = sprintf('./Example/S_outparAsymPar%d_%d_%d',...
        n_species,inet,isubnet);
    load(loadS);
    
    for i = 1:length(S_outpar)
        
        param = (isubnet-1)*length(S_outpar)+i;
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
    
    for j = 1:length(S_outpar)
        
        if  sol{j}.maxjackpot == 1
            if sol{j}.desc == 1
                if sol{j}.rightskewed == 1
                    if sol{j}.unimodal == 1
                        rare_par(count) = (isubnet-1)*10+j;
                        count = count + 1;
                    end
                end
            end
        end
    end
    
    %         save_sol = sprintf('/Volumes/MELANOMAII/Example/solAsymPar%d_%d_%d',n_species,inet,isubnet);
    save_sol = sprintf('./Example/solAsymPar%d_%d_%d',n_species,inet,isubnet);
    save(save_sol,'sol');
    
end

if exist('rare_par')
    %         saverarepar = sprintf('/Volumes/MELANOMAII/Example/rare_parAsymPar%d_%d',n_species,inet);
    saverarepar = sprintf('./Example/rare_parAsymPar%d_%d',n_species,inet);
    save(saverarepar,'rare_par')
    rare_par
else
    rare_par = [];
    %         saverarepar = sprintf('/Volumes/MELANOMAII/Example/rare_parAsymPar%d_%d',n_species,inet);
    saverarepar = sprintf('./Example/rare_parAsymPar%d_%d',n_species,inet);
    save(saverarepar,'rare_par')
    'No rare_par exists'
end

%% How many simulations do not show coordinatd high expression/unimodal right-skewed expression distributions

n_species = 2;
Net = 2;
Subnet = 1;

for inet = 1:Net
    inet
    count21{inet} = 0;
    count23{inet} = 0;
    for isubnet = 1:Subnet
        %         loadsol = sprintf('/Volumes/MELANOMAII/Example/sol%d_%d_%d',...
        %             n_species,inet,isubnet);
        loadsol = sprintf('./Example/sol%d_%d_%d',...
            n_species,inet,isubnet);
        load(loadsol)
        
        for j = 1:length(sol)
            if sol{j}.maxjackpot == 0
                count21{inet} = count21{inet} + 1;
            end
            if sol{j}.rightskewed == 0 || sol{j}.unimodal == 0
                count23{inet} = count23{inet} +1;
            end
        end
    end
end

% save('/Volumes/MELANOMAII/Example/count21','count21')
% save('/Volumes/MELANOMAII/Example/count23','count23')
save('./Example/count21','count21')
save('./Example/count23','count23')

