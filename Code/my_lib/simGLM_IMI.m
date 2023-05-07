function [tsp,sps,Itot,Istm] = simGLM_IMI(glmprs, Stim, rndseed)
rng(rndseed);
glmprs.nlfun = @exp;

% --------------- Check Inputs ----------------------------------
nbinsPerEval = 60;  % Default number of bins to update for each spike
% dt = glmprs(1).dtSp; % bin size for simulation

slen = size(Stim,1); % length of stimulus
rlen = slen;  % length of binned spike response

Istm = glmprs.kdc;
Itot = Istm; % total filter output
    
% --------------- Set up simulation dynamics variables ---------------
nsp = 0; % number of spikes
sps = zeros(rlen,1); % sparse spike time matrix
jbin = 1; % current time bin
tspnext = exprnd(1);  % time of next spike (in rescaled time)
rprev = 0;  % Integrated rescaled time up to current point
lastspk = 0;
last2spk = 0;
last3spk = 0;
last4spk = 0;

% --------------- Append zeros ---------------------------------------
if isfield(glmprs, 'ih1')
glmprs.ih1 = [glmprs.ih1; zeros(300, 1)];
end
if isfield(glmprs, 'ih2')
    glmprs.ih2 = [glmprs.ih2; zeros(300, 1)];
end
if isfield(glmprs, 'ih3')
    glmprs.ih3 = [glmprs.ih3; zeros(300, 1)];
end
if isfield(glmprs, 'ih4')
    glmprs.ih4 = [glmprs.ih4; zeros(300, 1)];
end
if isfield(glmprs, 'ih1h1h1')
    glmprs.ih1h1h1 = [glmprs.ih1h1h1; zeros(800, 1)];
end
% --------------- Run dynamics ---------------------------------------
while jbin <= rlen
    iinxt = jbin:min(jbin+nbinsPerEval-1,rlen);
    if isfield(glmprs, 'ih4')
        rrnxt = glmprs.nlfun(Itot(iinxt) + glmprs.ih1(iinxt - lastspk)...
             + glmprs.ih2(iinxt - last2spk) + glmprs.ih3(iinxt - last3spk) + glmprs.ih4(iinxt - last4spk) );    
    elseif isfield(glmprs, 'ih3')
        rrnxt = glmprs.nlfun(Itot(iinxt) + ...
            glmprs.ih1(iinxt - lastspk) + glmprs.ih2(iinxt - last2spk) + glmprs.ih3(iinxt - last3spk) ); 
    elseif isfield(glmprs, 'ih2')
        rrnxt = glmprs.nlfun(Itot(iinxt) + glmprs.ih1(iinxt - lastspk) + glmprs.ih2(iinxt - last2spk)); 
    elseif isfield(glmprs, 'ih1')
        rrnxt = glmprs.nlfun(Itot(iinxt) + glmprs.ih1(iinxt - lastspk)); 
    elseif isfield(glmprs, 'ih1h1')
        rrnxt = glmprs.nlfun(Itot(iinxt) + glmprs.ih12(iinxt - lastspk) + glmprs.ih12(iinxt - last2spk)); 
    elseif isfield(glmprs, 'ih1h1h1')
        rrnxt = glmprs.nlfun(Itot(iinxt) + glmprs.ih1h1h1(iinxt - lastspk) + glmprs.ih1h1h1(iinxt - last2spk) ...
            + glmprs.ih1h1h1(iinxt - last3spk));       
    else
        disp('does not contain history dpendency function')
    end
    
    rrcum = cumsum(rrnxt)+rprev; % integrated cond intensity
    
    % No spike in this window
    if (tspnext >= rrcum(end)) 
        jbin = iinxt(end)+1;
        rprev = rrcum(end);
        
    % Spike! ! !   
    else   
        ispk = iinxt(find(rrcum>=tspnext, 1, 'first')); % time bin where spike occurred
        nsp = nsp+1;
        sps(ispk) = 1; % spike time
        tspnext = exprnd(1);  % draw next spike time
        rprev = 0; % reset integrated intensity
        jbin = ispk+1;  % Move to next bin
                        
        %---- update synthesised firing rate ----
        iiPostSpk = (lastspk+1):ispk;        
        if isfield(glmprs, 'ih3')
            Itot(iiPostSpk) = Itot(iiPostSpk) + glmprs.ih1(iiPostSpk - lastspk)... 
                + glmprs.ih2(iiPostSpk - last2spk) + glmprs.ih3(iiPostSpk - last3spk);
        else
            disp('does not contain history dpendency function')
        end
        % --  Update # of samples per iter ---
        muISI = jbin/nsp;
        nbinsPerEval = max(20, round(1.5*muISI)); 
        
        last4spk = last3spk;
        last3spk = last2spk;
        last2spk = lastspk;
        lastspk = ispk;
            
    end
end
tsp = find(sps>0);


























