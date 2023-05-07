function [tsp,sps,Itot,Istm, In, Out] = simGLM_INTENSITY_bias(glmprs, Stim, rndseed)
%% concatenate intensity functions
kdc = ones(length(Stim), 1)*glmprs(1).dc;
for cc = 1:length(glmprs)
kdc(glmprs(cc).valid_point_range) = glmprs(cc).kdc(glmprs(cc).valid_point_range);
end

% figure('Position', [300, 300, 1400, 500]);
% plot(kdc)
% hold on
%% start simulation
rng(rndseed);

nlfun = @exp;
% T_wid = glmprs(1).T_wid;
% T_shift = glmprs(1).T_shift;
% T_wid_width = T_wid(end);
vv = {glmprs.valid_point_range};
% --------------- Check Inputs ----------------------------------
nbinsPerEval = 100;  % Default number of bins to update for each spike
dt = glmprs(1).dtSp; % bin size for simulation

slen = size(Stim,1); % length of stimulus
rlen = slen;  % length of binned spike response
hlen = size(glmprs(1).ih,1); % length of post-spike filter

% -------------  Compute filtered resp to signal ----------------
% Istm = glmprs(2).kdc;
Istm = kdc; 
Itot = Istm; % total filter output

% --------------- Set up simulation dynamics variables ---------------
nsp = 0; % number of spikes
sps = zeros(rlen,1); % sparse spike time matrix
jbin = 1; % current time bin
tspnext = exprnd(1);  % time of next spike (in rescaled time)
rprev = 0;  % Integrated rescaled time up to current point

In = [];
Out = [];
last_ispk = 0;

% --------------- Run dynamics ---------------------------------------
while jbin <= rlen
    iinxt = jbin:min(jbin+nbinsPerEval-1,rlen);
    rrnxt = nlfun(Itot(iinxt)); %----- *dt; % Cond Intensity
    rrcum = cumsum(rrnxt)+rprev; % integrated cond intensity
    if (tspnext >= rrcum(end)) % No spike in this window
        jbin = iinxt(end)+1;
        rprev = rrcum(end);
    else   % Spike!
        ispk = iinxt(find(rrcum>=tspnext, 1, 'first')); 
        % randomly select the boundary caused by the discrete bins.         
        tmpspk = find(rrcum>=tspnext, 1, 'first');
        sumright = rrcum(tmpspk);        
        tauout = sumright;
        
        if tmpspk > 1         
            sumleft = rrcum(tmpspk-1);
            ppleft = (sumright - tspnext)/(sumright - sumleft);
            
            if ppleft > 0.5
                ispk = ispk - 1;
                tauout = sumleft;
                if sumleft < 0.001
                    sumleft
                end
            end
        end
        
        nsp = nsp+1;
        sps(ispk) = 1; % spike time
        
        mxi = min(rlen, ispk+hlen); % max time affected by post-spike kernel
        iiPostSpk = ispk+1:mxi; % time bins affected by post-spike kernel
        if ~isempty(iiPostSpk)
        %-------  select which post-history filter should use -------
        idx = cellfun(@(x) nnz(ispk == x) > 0 ,vv);
        seg_id = find(idx);        
        if isempty(seg_id), seg_id = 1;end
        %----------------------------------------------------------
        Itot(iiPostSpk) = Itot(iiPostSpk)+glmprs(seg_id).ih(1:mxi-ispk);
        end
         
        tmp = sum(exp( Itot((last_ispk+1):ispk) )); 
        In = [In, tspnext];
        Out = [Out, tauout];
        last_ispk = ispk;
        
        tspnext = exprnd(1);  % draw next spike time
        rprev = 0; % reset integrated intensity
        jbin = ispk+1;  % Move to next bin
        % --  Update # of samples per iter ---
        muISI = jbin/nsp;
        nbinsPerEval = max(20, round(1.5*muISI)); 
    end
end
tsp = find(sps>0); %-----*dt;




















