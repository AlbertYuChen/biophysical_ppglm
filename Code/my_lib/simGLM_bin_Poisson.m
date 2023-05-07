function [tsp, sps, log_lambda, baseline] = simGLM_bin_Poisson(glmprs, Stim, rndseed)
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
hlen = size(glmprs.ih,1); % length of post-spike filter

% --------------- Set up simulation dynamics variables ---------------
sps = zeros(hlen+rlen,1); % append hlen zeros before the spike train. 
log_lambda = baseline;

for t = 1:rlen
    log_lambda(t) = dot( glmprs.ih, sps(t + hlen - (1:hlen)) ) + baseline(t);
    yy = poissrnd(exp( log_lambda(t) ));
    if yy ~= 0
        sps(t + hlen) = 1;
    end
end

sps = sps(hlen+1:end);
tsp = find(sps~=0);



















