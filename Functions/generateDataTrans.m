%This script produces Gillespie simulations according to random
%(latin hyper cube sampled) parameter sets for all weakly-connected
%non-isomorphic symmetric networks of a particular size (defined by
%n_species).
%
%INPUT:
%
%nruns:             number of different runs (= number of parameter sets)
%n_species:         number of network size
%maxgillespie:      number of Gillespie simulated time units
%gen:               'yes' = generate new parameters (or 'no' = load parameters)
%init:              initial values for all species and initial values bursts
%                   (where 0 = 'off' and 1 = 'on')
%data:              if gen = 'no', which parameter matrix to use
%type:              normal - normal parameter space; constrained - constrained
%                   parameter space; asym - if asymmetric architecture
%
%OUTPUT:
%
%S_outpar:           matrix of size 3*n_species times maxgillespie, where
%                   the columns give the gene produt count and DNA state
%                   (on/off) per time

%%
nruns = 1;
n_species = 5;
maxgillespie = 1000000;
gen = 'no';
init.Bon(1:n_species) = 0;
init.spec(1:n_species) = 20;
type = 'normal';

%overall specifications (for up to 10 species)
set_spec = cell(1,10);
set_P = cell(1,10);
set_Bon = cell(1,10);
set_Boff = cell(1,10);
set_prod = cell(1,10);
set_prodP = cell(1,10);
set_proddiff = cell(1,10);
set_deg = cell(1,10);
set_degP = cell(1,10);
set_onbasal = cell(1,10);
set_ondep = cell(1,10);
set_off = cell(1,10);

for iname = 1:10
    
    set_spec{iname} = sprintf('X%d', iname);   
    set_P{iname} = sprintf('P%d', iname); %species
    set_Bon{iname} = sprintf('B%d_on', iname);          %burst 'on'
    set_Boff{iname} = sprintf('B%d_off', iname);        %burst 'off'
    set_prod{iname} = sprintf('prod%d', iname);    
    set_prodP{iname} = sprintf('prodP%d', iname);       %production rate
    set_proddiff{iname} = sprintf('proddiff%d', iname); %difference in production rate between burst 'on' and burst 'off' (> 1)
    set_deg{iname} = sprintf('deg%d', iname);           %degradation rate
    set_degP{iname} = sprintf('degP%d', iname);     
    set_onbasal{iname} = sprintf('onbasal%d', iname);   %basal on-rate of burst
    set_ondep{iname} = sprintf('ondep%d', iname);       %additional on-rate due to dependency to other node
    set_off{iname} = sprintf('off%d', iname);           %off rate of burst
    
end

%set the initial values for B_off
for isetoff = 1: n_species
    if init.Bon(isetoff) == 1
        init.Boff(isetoff) = 0;
    else
        init.Boff(isetoff) = 1;
    end
end

%load all weakly-connected non-isomorphic symmetric networks of size
%n_species

% load_M_iso = sprintf('/Volumes/MELANOMA/Data/M_iso%d',n_species);
load_M_iso = sprintf('M_iso%d',n_species);
load(load_M_iso)

if isequal(gen,'yes') == 1                  %generate new parameters by latin hypercube sampling method
    rng('shuffle');                         %receive different values upon new start
    
    min_range = [repmat(0.01,1,1),...       %production rate
        repmat(0.001,1,1),...               %degradation rate
        repmat(0.001,1,1),...               %basal on-rate of burst
        repmat(0.1,1,1),...                 %Hill coefficient n (Hill function)
        repmat(0.1,1,1),...                 %additional on-rate due to dependency to other node
        repmat(0.01,1,1),...                %off rate of burst
        repmat(2,1,1)];                     %proddiff
    
    max_range = [repmat(1,1,1),...          %production rate
        repmat(0.1,1,1),...                 %degradation rate
        repmat(0.1,1,1),...                 %basal on-rate of burst
        repmat(10,1,1),...                  %Hill coefficient n (Hill function)
        repmat(1,1,1),...                   %additional on-rate due to dependency to other node
        repmat(0.1,1,1),...                 %off rate of burst
        repmat(100,1,1)];                   %proddiff
    
    latinhyp = lhsdesign_modified(nruns, min_range, max_range);
    
    latinhyp_prod = repmat(latinhyp(:,1),1,n_species);
    latinhyp_deg = repmat(latinhyp(:,2),1,n_species);
    latinhyp_onbasal = repmat(latinhyp(:,3),1,n_species);
    latinhyp_n = repmat(latinhyp(:,4),1,n_species);
    latinhyp_ondep = repmat(latinhyp(:,5),1,n_species);
    latinhyp_off = repmat(latinhyp(:,6),1,n_species);
    latinhyp_proddiff = repmat(latinhyp(:,7),1,n_species);
    
    k = repmat(0.95*latinhyp(:,1)./latinhyp(:,2).*latinhyp(:,7),1,n_species); %0.95 of stationary on-state of system
    
else
%     load('/Volumes/MELANOMA/Data/Data1000')  %load parameter set
    load('Data1000')  %load parameter set
    
    R_outpar_par = Data1000(968,:);

    latinhyp_prod = repmat(R_outpar_par(:,1),1,n_species);
    latinhyp_deg = repmat(R_outpar_par(:,2),1,n_species);
    latinhyp_onbasal = repmat(R_outpar_par(:,3),1,n_species);
    latinhyp_n = repmat(R_outpar_par(:,4),1,n_species);
    latinhyp_ondep = repmat(R_outpar_par(:,5),1,n_species);
    latinhyp_off = repmat(R_outpar_par(:,6),1,n_species);
    latinhyp_proddiff = repmat(R_outpar_par(:,8),1,n_species);
    a = 10;
    b = 10;
    k = repmat(0.95*latinhyp_prod(1)^2/latinhyp_deg(1)^2*latinhyp_proddiff(1)*a/b,1,n_species);
end


for istruc = 3
    
    fclose all;
    clear functions
    
    istruc
    
    clear T_outpar S_outpar
    
    T_outpar = cell(1,nruns);
    S_outpar = cell(1,nruns);
    S = zeros(nruns,4*n_species);
    R = zeros(nruns,10*n_species);
    P = zeros(nruns,6*n_species);
    
    for iruns = 1:nruns
        
        fclose all;
        clear functions
        clearvars -except nruns maxgillespie n_species init gen...
            set_spec set_P set_Bon set_Boff set_prod set_prodP set_proddiff set_deg set_degP...
            set_onbasal set_ondep set_off M_iso nstruc...
            latinhyp_prod latinhyp_deg latinhyp_onbasal latinhyp_n...
            latinhyp_ondep latinhyp_off latinhyp_proddiff k...
            iruns istruc S P R type a b
        
        parmat = M_iso{istruc};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %create txt file for specific network/sampled parameter set
        %         fileID = fopen('/Volumes/MELANOMAII/Example/gillespie_bursts.txt','w');
%         fileID = fopen('/Volumes/MELANOMAII/Revisions/gillespie_bursts.txt','w');
        fileID = fopen('gillespie_bursts.txt','w');
        
        fprintf(fileID,'%d\n', maxgillespie);
        
        for prod_txt = 1:n_species
            fprintf(fileID,'%s = %f : = %s \n', set_prod{prod_txt}, latinhyp_prod(iruns,prod_txt), set_spec{prod_txt}) ;
        end
        
        for prod_txt = 1:n_species
            fprintf(fileID,'%s = %f : %s = %s + %s \n', set_prodP{prod_txt}, a*latinhyp_prod(iruns,prod_txt), set_spec{prod_txt}, set_spec{prod_txt}, set_P{prod_txt}) ;
        end
        
        for deg_txt = 1:n_species
            fprintf(fileID,'%s = %f : %s = \n', set_deg{deg_txt}, latinhyp_deg(iruns,deg_txt), set_spec{deg_txt}) ;
        end
        
        for deg_txt = 1:n_species
            fprintf(fileID,'%s = %f : %s = \n', set_degP{deg_txt}, b*latinhyp_deg(iruns,deg_txt), set_P{deg_txt}) ;
        end
        
        for ondep_txt = 1:n_species
            fprintf(fileID,'%s = %f : %s = %s\n', set_ondep{ondep_txt}, latinhyp_ondep(iruns,ondep_txt), set_Boff{ondep_txt}, set_Bon{ondep_txt});
        end
        
        for on_txt = 1:n_species
            fprintf(fileID,'%s = %f : %s = %s\n', set_off{on_txt}, latinhyp_off(iruns,on_txt), set_Bon{on_txt}, set_Boff{on_txt});
        end
        
        fprintf(fileID,'%s\n', 'Production difference');
        
        for proddiff_txt = 1:n_species
            fprintf(fileID,'%s = %d\n', set_proddiff {proddiff_txt}, latinhyp_proddiff(iruns,proddiff_txt));
        end
        
        fprintf(fileID,'%s\n', 'Basal values');
        
        for onbasak_txt = 1:n_species
            fprintf(fileID,'%s = %d\n', set_onbasal{onbasak_txt}, latinhyp_onbasal(iruns,onbasak_txt));
        end
        
        fprintf(fileID,'%s\n', 'Hill function k');
        
        for spec_txt = 1:n_species
            fprintf(fileID,'k%s = %d\n', set_spec{spec_txt}, k(iruns,spec_txt));
        end
        
        fprintf(fileID,'%s\n', 'Hill function n');
        
        for n_txt = 1:n_species
            fprintf(fileID,'n%s = %d\n', set_spec{n_txt}, latinhyp_n(iruns,n_txt));
        end
        
        fprintf(fileID,'%s\n', 'Initial values');
        
        for initspec_txt = 1:n_species
            fprintf(fileID,'%s = %d\n', set_spec{initspec_txt}, init.spec(initspec_txt));
        end
        
        for initspec_txt = 1:n_species
            fprintf(fileID,'%s = %d\n', set_P{initspec_txt}, init.spec(initspec_txt));
        end
        
        for initB_txt = 1:n_species
            fprintf(fileID,'%s = %d\n', set_Boff{initB_txt}, init.Boff(initB_txt));
            fprintf(fileID,'%s = %d\n', set_Bon{initB_txt}, init.Bon(initB_txt));
        end
        
        fprintf(fileID,'%s\n', 'Network');
        for net_txt = 1:n_species
            fprintf(fileID,'%d\n',parmat(net_txt,:));
        end
        fclose(fileID);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fclose all;
        clear functions
        
        %         tok = make_param_mex_bursts('/Volumes/MELANOMAII/Example/gillespie_bursts.txt',...
        %             '/Volumes/MELANOMAII/Example/gillespie_bursts',n_species,iruns);

        tok = make_param_mex_bursts_trans('gillespie_bursts.txt',...
            'gillespie_bursts',n_species,iruns);
        
        gillespie_burstsparams  %load parameters set by txt
        
        S(iruns,:) = transpose(species);
        R(iruns,:) = transpose(rates);
        P(iruns,:) = transpose(propensity);
        
    end
    
    for kmem = 1:1
        
        for jruns = 1:1
            %         for jruns = 1:nruns
            
            par_spec = S(jruns,:);
            par_rates = R(jruns,:);
            par_prop = P(jruns,:);
            
            [times,savespecies] = gillespie_burstshistomex(0,par_spec,par_rates,par_prop,sum(clock*100),maxgillespie,maxgillespie);
            
            S_outpar{jruns} = savespecies;
            
            times = [];
            savespecies = [];
            
        end
        
        S_save = sprintf('S_outparTransFaster%d_%d_%d',n_species,istruc,kmem);
        
        save(S_save,'S_outpar','-v7.3');
        
    end
end

% %% check if rare behavior 
% 
% 
% % load('/Volumes/MELANOMAII/Revisions/S_outpar1000Points3_2_1')
% 
% %%%%%%%%FILL OUT%%%%%%%%%%%
% n_species = 3;
% Subnet = 1;
% inet = 2;
% doplot = 'no';
% tcell = 100000;                          %time per 'cell'
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % load('/Volumes/MELANOMAII/Example/Data50')
% load('/Volumes/MELANOMA/Data/Data1000')
% thres = Data1000(26,1)./Data1000(26,2).*Data1000(26,8)*0.8;
% 
% for isubnet = 1:Subnet
%     
%     clear sol
%     
%     for i = 1:length(S_outpar)
%         
%         [maxjackpot_sol,desc_sol,rightskewed_sol,unimodal_sol,samp_sol,samp_time_sol,rand_time_sol] ...
%             = analyzeQual_revision(n_species,tcell,thres,S_outpar{i});
%         sol{i}.maxjackpot = maxjackpot_sol;
%         sol{i}.desc = desc_sol;
%         sol{i}.rightskewed = rightskewed_sol;
%         sol{i}.unimodal = unimodal_sol;
%         sol{i}.samp = samp_sol;
%         sol{i}.samp_time = samp_time_sol;
%         sol{i}.time = rand_time_sol;
%     end
%     
% end
