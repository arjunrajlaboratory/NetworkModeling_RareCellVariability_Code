%function to produce Latin hypercube sampled parameter sets in all seven
%independent parameters (r_prod, r_deg, r_on, n, r_add, r_off, d) and
%create the dependent parameter k (where k = 0.95*r_prod/r_deg*d)
%
%INPUT:
%
%nruns:         number of parameter sets sampled
%n_species:     number of parameter multiples (eg. if 3 node network and
%               one requires a different parameter set per node: n_species = 3; 
%               
%type:          normal - normal parameter space; constrained - constrained
%               parameter space
%
%OUTPUT:
%
%D:             nruns * n_species*8 - matrix of nruns different LHS sampled
%               parameter sets for n_species

%%
function D = LHSsamplingAsymPar(nruns, n_species)

rng('shuffle');                             %receive different values upon new start

min_range = [repmat(0.01,1,n_species),...           %production rate
    repmat(0.001,1,n_species),...               %degradation rate
    repmat(0.001,1,n_species),...               %basal on-rate of burst
    repmat(5,1,n_species),...                 %Hill coefficient n (Hill function)
    repmat(0.2,1,n_species),...                 %additional on-rate due to dependency to other node
    repmat(0.01,1,n_species),...                %off rate of burst
    repmat(2,1,n_species)];                     %proddiff

max_range = [repmat(1,1,n_species),...      %production rate
    repmat(0.1,1,n_species),...             %degradation rate
    repmat(0.1,1,n_species),...             %basal on-rate of burst
    repmat(10,1,n_species),...              %Hill coefficient n (Hill function)
    repmat(0.4,1,n_species),...               %additional on-rate due to dependency to other node
    repmat(0.1,1,n_species),...             %off rate of burst
    repmat(100,1,n_species)];               %proddiff

latinhyp = lhsdesign_modified(nruns,min_range,max_range);

latinhyp_prod = latinhyp(:,1:n_species);
latinhyp_deg = latinhyp(:,n_species+1:2*n_species);
latinhyp_onbasal = latinhyp(:,2*n_species+1:3*n_species);
latinhyp_n = latinhyp(:,3*n_species+1:4*n_species);
latinhyp_ondep = latinhyp(:,4*n_species+1:5*n_species);
latinhyp_off = latinhyp(:,5*n_species+1:6*n_species);
latinhyp_proddiff = latinhyp(:,6*n_species+1:7*n_species);

k = 0.95*latinhyp(:,1:n_species)./latinhyp(:,n_species+1:2*n_species).*...
    latinhyp(:,6*n_species+1:7*n_species); %0.95 of stationary on-state of system

D = [latinhyp_prod, latinhyp_deg, latinhyp_ondep, latinhyp_off, latinhyp_proddiff,...
    latinhyp_onbasal, k, latinhyp_n];

end