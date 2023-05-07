function [tsp,sps,Itot,Istm] = simGLM_PosNeg(glmprs, Stim, StimExt, StimExt2, rndseed)
rng(rndseed);

% stim_smooth = Gaussian_Filter_YC(Stim, 4);
% stim_diff = [0; diff(stim_smooth)]*8;
% stim_diff_neg = stim_diff;
% stim_diff_neg(stim_diff_neg >0) = 0;

% stim_neg = Gaussian_Filter_YC(Stim, 6);
% stim_neg(stim_neg>=1) = 1;
% stim_neg = stim_neg - 1;
% stim_neg_sq = -stim_neg.^2;
%----------------------------------
glmprs.nlfun = @exp;
upSampFactor = glmprs.dtStim/glmprs.dtSp; % number of spike bins per Stim bin
assert(mod(upSampFactor,1) == 0, 'dtStim / dtSp must be an integer');

% --------------- Check Inputs ----------------------------------
nbinsPerEval = 100;  % Default number of bins to update for each spike

slen = size(Stim,1); % length of stimulus
rlen = slen*upSampFactor;  % length of binned spike response
hlen = size(glmprs.ih,1); % length of post-spike filter

% -------------  Compute filtered resp to signal ----------------
% IstmPos = sameconv(stim_pos, glmprs.posk);  
% IstmNeg = sameconv(stim_neg, glmprs.negk);  
% Istm = IstmPos + IstmNeg;

IstmPos = sameconv(Stim, glmprs.k);  
IstmExt1 = sameconv(StimExt, glmprs.Extk);  
IstmExt2 = sameconv(StimExt2, glmprs.Extk2);  
Istm = IstmPos + IstmExt1 + IstmExt2;
Istm = Istm + glmprs.dc; % upsample stim current
Itot = Istm; 

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


























