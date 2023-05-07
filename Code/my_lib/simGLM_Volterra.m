function [tsp,sps,Itot,Istm] = simGLM_PosNeg(glmprs, Stim, StimNln)
% [tsp,sps,Itot,Ispk] = simGLM(glmprs,Stim)
% 
% Compute response of glm to stimulus Stim.
%
% Uses time rescaling instead of Bernouli approximation to conditionally
% Poisson process
%
% Dynamics:  Filters the Stimulus with glmprs.k, passes this through a
% nonlinearity to obtain the point-process conditional intensity.  Add a
% post-spike current to the linear input after every spike.
%
% Input: 
%   glmprs - struct with GLM params, has fields 'k', 'h','dc' for params
%              and 'dtStim', 'dtSp' for discrete time bin size for stimulus
%              and spike train (in s).
%     Stim - stimulus matrix, with time running vertically and each
%              column corresponding to a different pixel / regressor.
% Output:
%   tsp - list of spike times (in s)
%   sps - binary matrix with spike times (at resolution dtSp).
%  Itot - summed filter outputs 
%  Istm - just the spike-history filter output

glmprs.nlfun = @exp;
upSampFactor = glmprs.dtStim/glmprs.dtSp; % number of spike bins per Stim bin
assert(mod(upSampFactor,1) == 0, 'dtStim / dtSp must be an integer');

% --------------- Check Inputs ----------------------------------
nbinsPerEval = 100;  % Default number of bins to update for each spike
dt = glmprs.dtSp; % bin size for simulation

slen = size(Stim,1); % length of stimulus
rlen = slen*upSampFactor;  % length of binned spike response
hlen = size(glmprs.ih,1); % length of post-spike filter

% -------------  Compute filtered resp to signal ----------------
DStim = [diff(Stim); 0]/glmprs.dtStim;
Istm = sameconv(Stim, glmprs.k);  % filter stimulus with k
IDstm = sameconv(DStim, glmprs.dk);  % filter stimulus with k
IstmNln = StimNln * glmprs.nln;
Inststm = Istm + glmprs.dc + IDstm + IstmNln; % upsample stim current
Itot = Inststm; 
    
% --------------- Set up simulation dynamics variables ---------------
nsp = 0; % number of spikes
sps = zeros(rlen,1); % sparse spike time matrix
jbin = 1; % current time bin
tspnext = exprnd(1);  % time of next spike (in rescaled time)
rprev = 0;  % Integrated rescaled time up to current point

% --------------- Run dynamics ---------------------------------------
while jbin <= rlen
    iinxt = jbin:min(jbin+nbinsPerEval-1,rlen);
    rrnxt = glmprs.nlfun(Itot(iinxt)); %--- *dt; % Cond Intensity
    rrcum = cumsum(rrnxt)+rprev; % integrated cond intensity
    if (tspnext >= rrcum(end)) % No spike in this window
        jbin = iinxt(end)+1;
        rprev = rrcum(end);
    else   % Spike!
        ispk = iinxt(find(rrcum>=tspnext, 1, 'first')); % time bin where spike occurred
        nsp = nsp+1;
        sps(ispk) = 1; % spike time
        mxi = min(rlen, ispk+hlen); % max time affected by post-spike kernel
        iiPostSpk = ispk+1:mxi; % time bins affected by post-spike kernel
        if ~isempty(iiPostSpk)
            Itot(iiPostSpk) = Itot(iiPostSpk)+glmprs.ih(1:mxi-ispk);
        end
        tspnext = exprnd(1);  % draw next spike time
        rprev = 0; % reset integrated intensity
        jbin = ispk+1;  % Move to next bin
        % --  Update # of samples per iter ---
        muISI = jbin/nsp;
        nbinsPerEval = max(20, round(1.5*muISI)); 
    end
end
tsp = find(sps>0); %----- *dt;


























