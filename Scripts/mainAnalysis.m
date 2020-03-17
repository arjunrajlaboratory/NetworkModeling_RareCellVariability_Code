%example of main analysis - for 50 parameter sets

% if ~exist('/Volumes/MELANOMAII/Example', 'dir')
%     mkdir /Volumes/MELANOMAII/Example
% end
if ~exist('./Example', 'dir')
    mkdir ./Example
end
addpath(genpath(pwd))

% addpath(genpath('/Volumes/MELANOMAII'))

%% create weakly-connected non-isomorphic symmetric networks (for up to 10 nodes)

clearvars;
clc;

n_species = 2;

M = zeros(1,n_species);

for icount = 1:n_species
    ind= [];
    for iM = 1:size(M,1)
        if sum(M(iM,:)) == icount - 1
            ind = [ind,iM];
        end
    end
    for jind = ind
        for icol = 1:n_species
            M_in = M(jind,:);
            M_in(icol) = 1;
            M(length(M)+1,:) = M_in;
        end
    end
    M = unique(M,'rows');
end

for j = 1:size(M,1)
    M_full{j}(1,:) = M(j,:);
    for ishift = 1:n_species -1
        M_full{j}(ishift+1,:) = circshift(M(j,:),ishift);
    end
end

count = 1;
for ifull = 1:length(M_full)
    G = digraph(M_full{ifull});
    weak_bins = conncomp(G,'Type','weak');
    if all(weak_bins == 1)
        M_full_con{count} = M_full{ifull};
        count = count + 1;
    end
end


for ifullunique = 1:length(M_full_con)
    if ifullunique == 1
        M_iso{1} = M_full_con{1};
    else
        G1 = digraph(M_full_con{ifullunique});
        countiso = 0;
        for iiso = 1:length(M_iso)
            G2 = digraph(M_iso{iiso});
            if isisomorphic(G1, G2) == 1
                break;
            else
                countiso = countiso + 1;
            end
        end
        if countiso == length(M_iso)
            M_iso{length(M_iso)+1} = M_full_con{ifullunique};
        end
    end
end

% M_iso_save = sprintf('/Volumes/MELANOMAII/Example/M_iso%d', n_species);
M_iso_save = sprintf('./Example/M_iso%d', n_species);
save(M_iso_save,'M_iso')

%% Generate LHS sampled parameter sets

clearvars
clc

%%%%%%%%%FILL OUT%%%%%%%%%%%
num_param = 50;
n_species = 1;              %as symmetric (same parameter sets used for all nodes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sample 50 parameter sets
Data50 = LHSsampling(num_param,n_species,'normal');

% save('/Volumes/MELANOMAII/Example/Data50','Data50')
save('./Example/Data50','Data50')

%% Generate data

%%%%%%%%%FILL OUT%%%%%%%%%%%
nruns = 50;                             %corresponds to number of parameter sets to be run
n_species = 2;                          %number of nodes
maxgillespie = 1000000;                 %number of Gillespie simulated time units
gen = 'no';                             %yes = generate new parameters (or 'no' = load parameters)
init.spec = repmat(20,1,n_species);     %initial values for all species (at t=0)
init.Bon = zeros(1,n_species);          %initial values for the burst (where 0 = 'off' and 1 = 'on')
data = 'Data50';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

generateData(nruns,n_species,maxgillespie,gen,init,data,'normal')

%% Population level

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
    
    %     loadS = sprintf('/Volumes/MELANOMAII/Example/S_outpar%d_%d_%d',...
    %         n_species,inet,isubnet);
    loadS = sprintf('./Example/S_outpar%d_%d_%d',...
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
    
    %     save_sol = sprintf('/Volumes/MELANOMAII/Example/sol%d_%d_%d',n_species,inet,isubnet);
    save_sol = sprintf('./Example/sol%d_%d_%d',n_species,inet,isubnet);
    save(save_sol,'sol');
    
end

if exist('rare_par')
    %     saverarepar = sprintf('/Volumes/MELANOMAII/Example/rare_par%d_%d',n_species,inet);
    saverarepar = sprintf('./Example/rare_par%d_%d',n_species,inet);
    save(saverarepar,'rare_par')
    rare_par
else
    rare_par = [];
    %     saverarepar = sprintf('/Volumes/MELANOMAII/Example/rare_par%d_%d',n_species,inet);
    saverarepar = sprintf('./Example/rare_par%d_%d',n_species,inet);
    save(saverarepar,'rare_par')
    'No rare_par exists'
end

%% Classification

clearvars
clc

%0 - stable low
%1 - uncorrelated transient high
%2 - rare transient correlaetd high
%3 - stably high

Net = 2;
n_species = 2;
Subnet = 1;
nruns = 50;
class = NaN(50*Net,4);

% load('/Volumes/MELANOMAII/Example/Data50')
load('./Example/Data50')
thres = Data50(:,1)./Data50(:,2).*Data50(:,8)*0.8;

for inet = 1:Net
    
    for i = 1:Subnet
        %         loadsol = sprintf('/Volumes/MELANOMAII/Example/sol%d_%d_%d',...
        %             n_species,inet,i);
        loadsol = sprintf('./Example/sol%d_%d_%d',...
            n_species,inet,i);
        load(loadsol);
        
        %         loadS = sprintf('/Volumes/MELANOMAII/Example/S_outpar%d_%d_%d',...
        %             n_species,inet,i);
        loadS = sprintf('./Example/S_outpar%d_%d_%d',...
            n_species,inet,i);
        load(loadS);
        
        for j = 1:length(sol)
            if sol{j}.maxjackpot == 0
                if sum(sum(S_outpar{j}(1:n_species,:)> thres(length(sol)*(i-1)+j))) == 0
                    class(inet*nruns-nruns+i*length(sol)-length(sol)+j,:) =...
                        [n_species, inet, length(sol)*(i-1)+j, 0];
                    
                elseif sum(sum(S_outpar{j}(1:n_species,:)> thres(length(sol)*(i-1)+j))) > 0
                    class(inet*nruns-nruns+i*length(sol)-length(sol)+j,:) =...
                        [n_species, inet, length(sol)*(i-1)+j, 1];
                    
                end
            else
                if sum(sol{j}.maxjackpot + sol{j}.desc + sol{j}.rightskewed + sol{j}.unimodal) == 4
                    class(inet*nruns-nruns+i*length(sol)-length(sol)+j,:) =...
                        [n_species, inet, length(sol)*(i-1)+j, 2];
                    
                elseif sum(sol{j}.maxjackpot + sol{j}.desc + sol{j}.rightskewed + sol{j}.unimodal) < 4
                    class(inet*nruns-nruns+i*length(sol)-length(sol)+j,:) =...
                        [n_species, inet, length(sol)*(i-1)+j, 3];
                    
                end
            end
        end
    end
end

% save('/Volumes/MELANOMAII/Example/class2','class')
save('./Example/class2','class')

%% Generate LHS sampled parameter sets on constrained parameter space

clearvars
clc

%%%%%%%%%FILL OUT%%%%%%%%%%%
num_param = 50;
n_species = 1;              %as symmetric (same parameter sets used for all nodes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sample 50 parameter sets
Data50C = LHSsampling(num_param,n_species,'constrained');

% save('/Volumes/MELANOMAII/Example/Data50C','Data50C')
save('./Example/Data50C','Data50C')

%% Generate data on constrained parameter space

clearvars
clc

%%%%%%%%%FILL OUT%%%%%%%%%%%
nruns = 50;                           %corresponds to number of parameter sets to be run
n_species = 2;                          %number of nodes
maxgillespie = 1000000;                 %number of Gillespie simulated time units
gen = 'no';                             %yes = generate new parameters (or 'no' = load parameters)
init.spec = repmat(20,1,n_species);     %initial values for all species (at t=0)
init.Bon = zeros(1,n_species);          %initial values for the burst (where 0 = 'off' and 1 = 'on')
data = 'Data50C';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%specify in function to only run eg. network 3_2 (and later 5_3)
generateData(nruns,n_species,maxgillespie,gen,init,data,'constrained')


%% Population level for constrained parameter space

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

% load('/Volumes/MELANOMAII/Example/Data50C')
load('./Example/Data50C')
thres = Data50C(:,1)./Data50C(:,2).*Data50C(:,8)*0.8;

for isubnet = 1:Subnet
    
    clear sol
    
    %     loadS = sprintf('/Volumes/MELANOMAII/Example/S_outparC%d_%d_%d',...
    %         n_species,inet,isubnet);
    loadS = sprintf('./Example/S_outparC%d_%d_%d',...
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
    
    %     save_sol = sprintf('/Volumes/MELANOMAII/Example/solC%d_%d_%d',n_species,inet,isubnet);
    save_sol = sprintf('./Example/solC%d_%d_%d',n_species,inet,isubnet);
    save(save_sol,'sol');
    
end

if exist('rare_par')
    %     saverarepar =
    %     sprintf('/Volumes/MELANOMAII/Example/rare_parC%d_%d',n_species,inet);
    saverarepar = sprintf('./Example/rare_parC%d_%d',n_species,inet);
    save(saverarepar,'rare_par')
    rare_par
else
    rare_par = [];
    %     saverarepar = sprintf('/Volumes/MELANOMAII/Example/rare_parC%d_%d',n_species,inet);
    saverarepar = sprintf('./Example/rare_parC%d_%d',n_species,inet);
    save(saverarepar,'rare_par')
    'No rare_par exists'
end

%% Burst Analysis

%%%%%%%%%FILL OUT%%%%%%%%%%%
n_species = 2;
inet = 2;
startSubnet = 1;
endSubnet = 1;
Subnet = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_analyzeBurstsAll(n_species, inet, startSubnet, endSubnet, Subnet)
