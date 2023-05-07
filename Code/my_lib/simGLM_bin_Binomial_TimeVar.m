function [tsp, sps, eta, baseline] = simGLM_bin_Binomial_TimeVar(glmprs, Stim, rndseed)
%% concatenate intensity functions
kdc = glmprs.kdc;
kdc(glmprs.valid_point_range) = glmprs.kdc(glmprs.valid_point_range);
baseline = kdc; 

% figure('Position', [300, 300, 1400, 500]);
% plot(kdc)
% hold on
% 
% figure
% plot(glmprs.ih)

%% start simulation
rng(rndseed);

slen = size(Stim,1); % length of stimulus
rlen = slen;  % length of binned spike response
hlen = size(glmprs.h,1); % length of post-spike filter

lowerBnd = log(eps('double')); 
upperBnd = -lowerBnd;
sigmoid = @(eta) 1 ./ (1 + exp(-constrain(eta,lowerBnd,upperBnd)));

% --------------- Set up simulation dynamics variables ---------------
sps = zeros(hlen+rlen,1); % append hlen zeros before the spike train. 
eta = baseline;

for t = 1:rlen
    eta(t) = dot( glmprs.h(:,t), sps(t + hlen - (1:hlen)) ) + baseline(t);
    yy = binornd(1, sigmoid( eta(t) ));
    if yy ~= 0
        sps(t + hlen) = 1;
    end
end

sps = sps(hlen+1:end);
tsp = find(sps~=0);


%-----------------------------------------------------------------------
function x = constrain(x,lower,upper)
% Constrain between upper and lower limits, and do not ignore NaN
x(x<lower) = lower;
x(x>upper) = upper;
















