%function to receive population levels for given simulations and used to
%determine the class of each simulation 
%
%
%INPUT:
%
%n_species:     number of nodes in the network of given simulation (eg. 2, 3, 5, 8)
%tcell:         time per cell (here: 1000 time units per cell)
%threshold:     threshold at which the cell is called to highly express the
%               gene (here: calculated as the 0.8*steady state of high state
%S_outpar:      output of simulation file - number of gene product counts per
%               time unit per gene and state of DNA per gene per time unit (on/off)
%param:         the parameter set for which the population level should be
%               calculated for (here: [1,1000])
%inet:          determine for which network architecture of the defined node size
%               the population level should be caclulated for (eg. for network size 3, 
%               there are 4 different network architectures, hence param /in [1,4])
%
%OUTPUT:
%
%maxjackpot:    1 - at least once throughout the simulation the number of
%                   highly expressed genes exceeds half the network size 
%                   (eg. for a 5 node network at least once 3 genes are highly 
%                   expressed at the same time) 
%               0 - else: not once throughout the whole simulation does the
%                   number of highly expresssed genes exceed half the network
%                   size (eg. for a 5 node network there are never more than
%                   2 genes highly expressed at the same time)
%desc:          1 - the histogram of number of highly expressed genes per snapshot 
%                   per cell is monotonically descreasing, where we take 
%                   tcell time units of the whole simulation to be one cell 
%                   and randomly determine one snapshot time point of that
%                   intervall
%               0 - the histogram of number of highly expressed genes per snapshot 
%                   per cell is not monotonically descreasing 
%rightskewed:   1 - the gene expression distribution per snapshot 
%                   per cell is right skewed according to the MATLAB function 
%                   'skewness'
%               0 - the gene expression distribution per snapshot 
%                   per cell is not right skewed according to the MATLAB function 
%                   'skewness'
%unimodal:      1 - the gene expression distrbution is 'unimodal': the
%                   minimum of the first quarter of the gene expression 
%                   distribution per snapshot per cell is greater than the 
%                   minimum of the last quarter of the gene expression dsitribution
%               0 - the gene expression distrbution is not 'unimodal'
%samp:          summarizes the gene product counts at the randomly
%               determined time points per cell (gives gene expression
%               distribution)
%samp_time:     the number of highly expressed genes per snapshot per cell
%               (gives histogram of number of highly expressed genes)
%rand_time:     randomly determined timepoint at which 'snapshot' is taken
%               per cell (\in [1,tcell])

%%
function [maxjackpot,desc,rightskewed,unimodal,samp,samp_time,rand_time] = ...
    analyzeQual_revision(n_species,tcell,threshold,S_outpar)

if length(threshold) ~= n_species
   threshold = repmat(threshold,1,n_species);
end

rng('shuffle')
rand_time = datasample(1:tcell,1);

for j = 1:n_species
    cells(j,1:length(S_outpar)) = S_outpar(j,:) > threshold(j);
end
all_cells = sum(cells);

samp = S_outpar(1:n_species,rand_time:tcell:end);

samp_time = zeros(1,length(S_outpar(1:n_species,rand_time:tcell:end)));

for k = 1:n_species
    for jsamp = 1:length(samp)
        samp_time1(k,jsamp) = samp(k,jsamp) > threshold(k);
    end
end
samp_time = sum(samp_time1);

[val,ind] = sort(histcounts(samp_time,-0.5:1:n_species+1),'descend');

if mod(n_species,2) == 1
    if max(all_cells) >= ceil(n_species/2)
        maxjackpot = 1;
    else
        maxjackpot = 0;
    end
else
    if max(all_cells) >= ceil(n_species/2) + 1
        maxjackpot = 1;
    else
        maxjackpot = 0;
    end
end

if isequal(0:n_species,ind-1) == 1
    desc = 1;
else
    desc = 0;
end

for ispec = 1:n_species
    N = histcounts(samp(ispec,:));
    if skewness(N) > 0
         rightskewed_spec(ispec) = 1;
    else
        rightskewed_spec(ispec) = 0;
    end
    quarter = ceil(length(N)/4);
    if min(N(1:quarter)) > max(N(end-quarter:end))
        unimodal_spec(ispec) = 1;
    else
        unimodal_spec(ispec) = 0;
    end
end

rightskewed = min(rightskewed_spec);
unimodal = min(unimodal_spec);


end