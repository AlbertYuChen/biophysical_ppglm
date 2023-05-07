function [tsp,sps,Itot,Istm, In, Out] = simGLM_Pillow_bias(glmprs, Stim, rndseed)
rng(rndseed);

glmprs.nlfun = @exp;

% --------------- Check Inputs ----------------------------------
nbinsPerEval = 100;  % Default number of bins to update for each spike
dt = glmprs.dtSp; % bin size for simulation

slen = size(Stim,1); % length of stimulus
rlen = slen;  % length of binned spike response
hlen = size(glmprs.ih,1); % length of post-spike filter

% -------------  Compute filtered resp to signal ----------------
% DStim = [diff(Stim); 0]/glmprs.dtStim;
% Istm = sameconv(Stim, glmprs.k);  % filter stimulus with k
% IDstm = sameconv(DStim, glmprs.dk);  % filter stimulus with k
% Inststm = Istm + glmprs.dc ; % upsample stim current

Istm = glmprs.kdc; 
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
    rrnxt = glmprs.nlfun(Itot(iinxt)); %------- *dt; % Cond Intensity
    rrcum = cumsum(rrnxt)+rprev; % integrated cond intensity
    if (tspnext >= rrcum(end)) % No spike in this window
        jbin = iinxt(end)+1;
        rprev = rrcum(end);
    else   % Spike! 
        ispk = iinxt(find(rrcum>=tspnext, 1, 'first')); % time bin where spike occurred
        % randomly select the boundary caused by the discrete bins. 
        % Add by Yu Chen July 8th 2017       
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
        
        
        tmp = sum(exp( Itot((last_ispk+1):ispk) ));
        In = [In, tspnext];
        Out = [Out, tauout];
        last_ispk = ispk;
        
              
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
tsp = find(sps>0); %----------- *dt;


























