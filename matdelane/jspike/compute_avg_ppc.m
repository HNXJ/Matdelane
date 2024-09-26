function [w_ppc,uw_ppc,thresh_ppc] = compute_avg_ppc(n_spikes,ppc_all)

% compute weighted ppc
nNeurons = size(ppc_all,1);
w_ppc = nansum ( (ppc_all .* n_spikes) ./ repmat(nansum(n_spikes,1),[nNeurons 1]) ,1); 
%--to compute a SEM one would have to take a bootstrap here 

% unweighted average 
uw_ppc = nanmean(ppc_all,1); 

% threshold of 50 spikes
ppc_all_thresh = ppc_all; 
ppc_all_thresh(n_spikes<50) = NaN;
thresh_ppc = nanmean(ppc_all_thresh,1); 
