% Author: Yu Chen
% Date: Aug 14th 2018 @ Beijing

function [tsp, sps, eta, baseline] = simGLM_bin_Binomial_slidingH(glmprs, Stim, rndseed)
%% concatenate intensity functions
baseline = glmprs.kdc;

%% start simulation
rng(rndseed);

slen = size(Stim,1); % length of stimulus
rlen = slen;  % length of binned spike response
hlen = size(glmprs.iht,1); % length of post-spike filter

blk_N = glmprs.blk_N;
block_range = glmprs.block_range;

% --------------- Set up simulation dynamics variables ---------------
sps = zeros(hlen+rlen,1); % append hlen zeros before the spike train. 
eta = baseline;

for t = 1:rlen
    idx = cellfun(@(x)ismember(t,x)>0, block_range);
    if sum(idx)==0; idx=1; end
    eta(t) = dot( glmprs.h{idx}, sps(t + hlen - (1:hlen)) ) + baseline(t);
    yy = binornd(1, sigmoid( eta(t) ));
    if yy ~= 0
        sps(t + hlen) = 1;
    end
end

sps = sps(hlen+1:end);
tsp = find(sps~=0);





























