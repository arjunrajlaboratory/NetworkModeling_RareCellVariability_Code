%This script produces Gillespie simulations according to random
%(latin hyper cube sampled) parameter sets for all weakly-connected
%non-isomorphic symmetric networks of a particular size (defined by
%n_species).

function S_outpar = generateDataHierarchical(nruns,n_species,maxgillespie,gen,data)

init.spec = repmat(20,1,n_species);
init.Bon = zeros(1,n_species); 

%overall specifications (for up to 10 species)
set_spec = cell(1,10);
set_Bon = cell(1,10);
set_Boff = cell(1,10);
set_prod = cell(1,10);
set_proddiff = cell(1,10);
set_deg = cell(1,10);
set_onbasal = cell(1,10);
set_ondep = cell(1,10);
set_off = cell(1,10);

for iname = 1:10
    
    set_spec{iname} = sprintf('X%d', iname);            %species
    set_Bon{iname} = sprintf('B%d_on', iname);          %burst 'on'
    set_Boff{iname} = sprintf('B%d_off', iname);        %burst 'off'
    set_prod{iname} = sprintf('prod%d', iname);         %production rate
    set_proddiff{iname} = sprintf('proddiff%d', iname); %difference in production rate between burst 'on' and burst 'off' (> 1)
    set_deg{iname} = sprintf('deg%d', iname);           %degradation rate
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
M_iso{1} = repmat([1 0 0 0 0],5,1);

nstruc = length(M_iso);                     %number of all possible networks considered

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
    R_outpar_par = data;
    latinhyp_prod = repmat(R_outpar_par(:,1),1,n_species);
    latinhyp_deg = repmat(R_outpar_par(:,2),1,n_species);
    latinhyp_onbasal = repmat(R_outpar_par(:,3),1,n_species);
    latinhyp_n = repmat(R_outpar_par(:,4),1,n_species);
    latinhyp_ondep = repmat(R_outpar_par(:,5),1,n_species);
    latinhyp_off = repmat(R_outpar_par(:,6),1,n_species);
    latinhyp_proddiff = repmat(R_outpar_par(:,8),1,n_species);
    
    k = repmat(R_outpar_par(:,7),1,n_species);
end


for istruc = 1
    
    fclose all;
    clear functions
    
    istruc
    
    clear T_outpar S_outpar
    
    T_outpar = cell(1,nruns);
    S_outpar = cell(1,nruns);
    S = zeros(nruns,3*n_species);
    R = zeros(nruns,8*n_species);
    P = zeros(nruns,4*n_species);
    
    for iruns = 1:nruns
        
        fclose all;
        clear functions
        clearvars -except nruns maxgillespie n_species init gen...
            set_spec set_Bon set_Boff set_prod set_proddiff set_deg...
            set_onbasal set_ondep set_off M_iso nstruc...
            latinhyp_prod latinhyp_deg latinhyp_onbasal latinhyp_n...
            latinhyp_ondep latinhyp_off latinhyp_proddiff k...
            iruns istruc S P R
        
        parmat = M_iso{istruc};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %create txt file for specific network/sampled parameter set
        fileID = fopen('gillespie_bursts.txt','w');
        fprintf(fileID,'%d\n', maxgillespie);
        
        for prod_txt = 1:n_species
            fprintf(fileID,'%s = %f : = %s \n', set_prod{prod_txt}, latinhyp_prod(iruns,prod_txt), set_spec{prod_txt}) ;
        end
        
        for deg_txt = 1:n_species
            fprintf(fileID,'%s = %f : %s = \n', set_deg{deg_txt}, latinhyp_deg(iruns,deg_txt), set_spec{deg_txt}) ;
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
        
        tok = make_param_mex_bursts('gillespie_bursts.txt','gillespie_bursts',n_species,iruns);
        
        
        gillespie_burstsparams  %load parameters set by txt
        
        S(iruns,:) = transpose(species);
        R(iruns,:) = transpose(rates);
        P(iruns,:) = transpose(propensity);
        
    end
    
    R_save = sprintf('R_outpar%d_%d',n_species,istruc);
    S_save = sprintf('S_outpar%d_%d',n_species,istruc);
    P_save = sprintf('P_outpar%d_%d',n_species,istruc);
    
    save(R_save,'R','-v7.3');
    save(S_save,'S','-v7.3');
    save(P_save,'P','-v7.3');
    
    for kmem = 1:10
        
        parfor jruns = 1:100
            %         for jruns = 1:nruns
            
            par_spec = S(100*kmem-100+jruns,:);
            par_rates = R(100*kmem-100+jruns,:);
            par_prop = P(100*kmem-100+jruns,:);
            
            [times,savespecies] = gillespie_burstshistomex(0,par_spec,par_rates,par_prop,sum(clock*100),maxgillespie,maxgillespie);
            
            S_outpar{jruns} = savespecies;
            
            savespecies = [];
            
        end
        
        S_save = sprintf('S_outparHierarchical%d_%d_%d',n_species,istruc,kmem);
        
        save(S_save,'S_outpar','-v7.3');
        
        clear S_outpar
    end
end

end
