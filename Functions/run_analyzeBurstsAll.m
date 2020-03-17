%function analysing bursting 
%
%INPUT:
%
%n_species:             number od nodes in network
%inet:                  number of network architecture (eg. for a 3 node
%                       network there are 4 architectures)
%startSubnet:           due to limited internal storage the simulations are
%                       further devided into subnetworks (eg. the 1000 
%                       simulations for one 3 node archticture are saved in 
%                       10 different subnetwork files - each contianing
%                       100 simulations), startSubnet defines at which
%                       subnetwork the analysis starts
%endSubnet:             defines at which subnetwork the analysis ends
%Subnet:                number of all subnetworks (eg. for 3 nodes networks
%                       Subnet = 10)
%
%OUTPUT:
%
%solBurstsOnInit:       structure M containing fields
%                           BurstOnJInit: median duration of bursts at entry 
%                                         time point for each simulation showing 
%                                         rare coordinated high states 
%                           BurstOnJInit_Data: durationf of all bursts at 
%                                         entry time point for each 
%                                         simulation showing rare coordinated 
%                                         high states
%                           BurstOnJNInit: median duration of bursts in baseline
%                                          time region for each simulation showing 
%                                          rare coordinated high states 
%                           BurstOnJNInit_Data: durationf of all bursts in baseline
%                                          time region for each simulation showing 
%                                          rare coordinated high states
%                           StatSig: statistical significance according to
%                                    kstest2 (for BurstOnJInit_Data and
%                                    BurstOnJNInit_Data) where 1 -
%                                    signifcant difference at level 0.05,
%                                    0 - no significant difference at level 0.05
%                                    per simulation showing rare coordinated high states
%                           P: corresponding p values to kstest2 per simulation 
%                              showing rare coordinated high states
%                           ks2Stat: ks2stat statistic - see MATLAB
%                                    function kstest2 for more information 
%                           Par: all simulations/parameters showing rare 
%                               coordinated high states for that subnetwork
%
%solBurstsTerm:       structure M containing fields
%                           BurstOnJExit_Data: durationf of all bursts in exit 
%                                         time region for each simulation
%                                         showing rare coordinated 
%                                         high states
%                           BurstOnJNExit_Data: durationf of all bursts in 
%                                         high time region but not in exit 
%                                         time region for each simulation
%                                         showing rare coordinated 
%                                         high states
%                           StatSigOn: statistical significance according to
%                                    kstest2 (for BurstOnJExit_Data and
%                                    BurstOnJNExit_Data) where 1 -
%                                    signifcant difference at level 0.05,
%                                    0 - no significant difference at level 0.05
%                                    per simulation showing rare coordinated high states
%                           POn: corresponding p values to kstest2 per simulation 
%                              showing rare coordinated high states
%                           ks2StatOn: ks2stat statistic - see MATLAB
%                                    function kstest2 for more information 
%                           Par: all simulations/parameters showing rare 
%                               coordinated high states for that subnetwork
%
%solEnterJack:          structure sol with fields
%                            Data: the times interval between two consecutive 
%                                  nodes exceeding the threshold for all 
%                                  simulation showing rare coordinated 
%                                  high states and all high states
%                            H: test decision for lilliefors test (null
%                               hypothesis: Data is exponentially distributed)
%                            P: corresponding p value
%                            K: test statistic kstat - see MATLAB function 
%                               lilliefors for more information 
%                            C: critical vlaue for test
%solExitJack:          structure sol with fields
%                            Data: the times interval between two consecutive 
%                                  nodes falling below the threshold for all 
%                                  simulation showing rare coordinated 
%                                  high states and all high states
%                            H: test decision for lilliefors test (null
%                               hypothesis: Data is exponentially distributed)
%                            P: corresponding p value
%                            K: test statistic kstat - see MATLAB function 
%                               lilliefors for more information 
%                            C: critical vlaue for test

%%
function  run_analyzeBurstsAll(n_species, inet, startSubnet, endSubnet, Subnet)

clearvars -except n_species startSubnet endSubnet inet Subnet
clc;
doplot = 'no';
data = 'Data50';


% rare = sprintf('/Volumes/MELANOMAII/Example/rare_par%d_%d',...
%     n_species, inet);

rare = sprintf('./Example/rare_par%d_%d',...
    n_species, inet);
load(rare)

load(data);
thres = Data50(:,1)./Data50(:,2).*Data50(:,8)*0.8;

for isubnet = startSubnet:endSubnet
    
%     loadS = sprintf('/Volumes/MELANOMAII/Example/S_outpar%d_%d_%d',n_species, inet, isubnet);
    loadS = sprintf('./Example/S_outpar%d_%d_%d',n_species, inet, isubnet);
    load(loadS);

    r1 = rare_par(rare_par > isubnet*length(S_outpar)-length(S_outpar));
    r2 = rare_par(rare_par <= isubnet*length(S_outpar)-length(S_outpar)+50/Subnet);
    rare_par_isubnet = intersect(r1,r2);
    
    analyzeBurstsInit(n_species, S_outpar, rare_par_isubnet, inet, thres, isubnet);

    analyzeBurstsTerm(n_species, S_outpar, rare_par_isubnet, inet, thres, isubnet);

    analyzeEnterJack(n_species, S_outpar, rare_par_isubnet, doplot, inet, thres, isubnet);
    
    analyzeExitJack(n_species, S_outpar, rare_par_isubnet, doplot, inet, thres,isubnet);
end
