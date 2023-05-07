function [tsp,sps,Itot,Istm] = simGLM_Isolink(glmprs, Stim, rndseed)
%% concatenate intensity functions
kdc = glmprs.kdc;
kdc(glmprs.valid_point_range) = glmprs.kdc(glmprs.valid_point_range);
Istm = kdc; 

% figure('Position', [300, 300, 1400, 500]);
% plot(kdc)
% hold on
% 
% figure
% plot(glmprs.ih)

%% start simulation
rng(rndseed);

dt = glmprs.dtSp; % bin size for simulation
slen = size(Stim,1); % length of stimulus
rlen = slen;  % length of binned spike response
hlen = size(glmprs.ih,1); % length of post-spike filter

% --------------- Set up simulation dynamics variables ---------------
sps = zeros(hlen+rlen,1); % append hlen zeros before the spike train. 
Itot = Istm;

% setup link function
last_zero_ind = find(glmprs.link_func(:,2) == 0, 1, 'last');
glmprs.link_func(1:(last_zero_ind-2),:) = []; % delete redundant zeros on the left of the link functions
glmprs.link_func(1,1) = glmprs.link_func(2,1)-1;
glmprs.link_func = [glmprs.link_func; [glmprs.link_func(end,1)+1, 1]]; % set the right bound to a larger value

x_range = round(glmprs.link_func(2,1), 3):dt:round(glmprs.link_func(end-1,1), 3);
index_search_link_func = interp1q(glmprs.link_func(:,1), glmprs.link_func(:,2), x_range');

left_bound = x_range(1);
right_bound = x_range(end);

%----------- test link function here -----------
% xin = (-8:0.001:8)';
% x = xin;
% x(x<left_bound) = left_bound;
% x(x>right_bound) = right_bound;
% search_index = round((x-left_bound)*1000)+1;
% link_out = index_search_link_func(search_index);
% 
% figure
% plot(xin, link_out)
%----------- test link function here -----------

for t = 1:rlen
    Itot(t) = dot( glmprs.ih, sps(t + hlen - (1:hlen)) ) + Istm(t);
    xin = Itot(t);
    % input to link function 
    x = xin;
    x(x<left_bound) = left_bound;
    x(x>right_bound) = right_bound;
    search_index = round((x-left_bound)*1000)+1;
    lambda = index_search_link_func(search_index);

    yy = poissrnd(lambda);
    if yy ~= 0
        sps(t + hlen) = 1;
    end
end

sps = sps(hlen+1:end);
tsp = find(sps~=0);



















