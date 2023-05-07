
function [link,dlink,ilink] = GLM_LINK_YC
dataClass = 'double';
tiny = realmin(dataClass)^.25; % keep fourth powers from under/overflowing

%% ---- Poisson link function -------
link0 = @(mu) log(mu);
dlink0 = @(mu) 1 ./ mu;
% keep mu = ilink(eta) in [tiny, 1/tiny];
lowerBnd = log(tiny); upperBnd = -lowerBnd;
ilink0 = @(eta) exp(constrain_YC(eta,lowerBnd,upperBnd));

%% ---- Bernoulli -------
link1 = @(mu) log(mu ./ (1-mu));
dlink1 = @(mu) 1 ./ (mu .* (1-mu));
% keep mu = ilink(eta) in (approx) [eps, (1-eps)];
lowerBnd = log(eps(dataClass)); upperBnd = -lowerBnd;
ilink1 = @(eta) 1 ./ (1 + exp(-constrain(eta,lowerBnd,upperBnd)));


%% ---- mu = exp(CC * x) ------- this is used to prove CC is not changing
% the performance of GLM framework, but the amplitude of the beta values
% CC = 0.2;
% link = @(mu) log(mu)/CC;
% dlink = @(mu) 1./ mu / CC;
% lowerBnd = log(tiny); upperBnd = -lowerBnd;
% ilink = @(eta)exp(CC*constrain_YC(eta,lowerBnd,upperBnd));


%% ---- mu = exp(x^2) -------
% link = @(mu) sqrt(log(mu));
% dlink = @(mu) 1 ./ mu ./ link(mu);
% % keep mu = ilink(eta) in (approx) [eps, (1-eps)];
% lowerBnd = log(eps(dataClass)); upperBnd = -lowerBnd;
% ilink = @(eta) exp(constrain(eta,lowerBnd,upperBnd).^2);

%% self defined
link = link0; 
dlink = dlink0; 
ilink = ilink0;

% link = @(mu) linkExt(mu); 
% dlink = @(mu) dlinkExt(mu); 
% lowerBnd = log(tiny); upperBnd = -lowerBnd;
% ilink = @(eta) ilinkExt(constrain(eta,lowerBnd,upperBnd));

% Coe = [-0.1453   -0.5648    1.1034   -0.1700];
% link = @(mu) linkExt(mu, Coe); 
% dlink = @(mu) dlinkExt(mu, Coe); 
% lowerBnd = log(tiny); upperBnd = -lowerBnd;
% ilink = @(eta) ilinkExt(constrain(eta,lowerBnd,upperBnd), Coe);


% t=[0.01:0.01:3]';
% figure
% plot(t, link(t))
% hold on
% plot(t, link0(t))
% plot(t, link1(t))
% grid on
% 
% figure
% plot(t, dlink(t))
% hold on
% plot(t, dlink0(t))
% plot(t, dlink1(t))
% grid on
% 
% t = [-3:0.01:2]';
% figure
% plot(t, ilink(t))
% hold on
% plot(t, ilink0(t))
% plot(t, ilink1(t))
% grid on


%% exp + straight line 
% function out = linkExt(mu)
% CC = 0;
% out = ones(size(mu));
% Ind1 = mu > exp(CC);
% Ind2 = mu <= exp(CC);
% 
% out(Ind1) = (mu(Ind1) - exp(CC))/exp(CC) + CC;
% out(Ind2) = log(mu(Ind2));
% 
% 
% function out = dlinkExt(mu)
% CC = 0;
% out = ones(size(mu));
% Ind1 = mu > exp(CC);
% Ind2 = mu <= exp(CC);
% 
% out(Ind1) = 1 / exp(CC);
% out(Ind2) = 1 ./ mu(Ind2);
% 
% 
% function out = ilinkExt(eta)
% CC = 0;
% out = ones(size(eta));
% Ind1 = eta > CC;
% Ind2 = eta <= CC;
% 
% out(Ind1) = (eta(Ind1)-CC)*exp(CC) + exp(CC);
% out(Ind2) = exp(eta(Ind2));

%% exp + straight line 
% function out = linkExt(mu)
% out = mu;
% 
% function out = dlinkExt(mu)
% out = mu;
% 
% Ind1 = mu > 0;
% Ind2 = mu <= 0;
% 
% out(Ind1) = ones(sum(Ind1),1);
% out(Ind2) = 9999999999999*ones(sum(Ind2),1);
% 
% 
% function out = ilinkExt(eta)
% out = eta;
% 
% Ind1 = eta > 0;
% Ind2 = eta <= 0;
% 
% out(Ind1) = eta(Ind1);
% out(Ind2) = 0;


%% exp shifted
% function out = linkExt(mu)
% SS = 0;
% out = log(mu+SS);
% 
% 
% function out = dlinkExt(mu)
% SS = 0;
% out = 1 ./ (mu+SS);
% 
% 
% function out = ilinkExt(eta)
% SS = 0;
% out = exp(eta) - SS;
% out(out<0) = 0;


%% fast decay + normal exp
function out = linkExt(mu, Coe)
CC = -1.5;
out = ones(size(mu));
Ind1 = mu <= exp(CC);
Ind2 = mu > exp(CC);

% out = log(mu);
% out(Ind1) = (mu(Ind1) - exp(CC))/exp(CC) + CC;
% out = ( log( (log(mu)-Coe(4))/Coe(1) ) - Coe(3)) /Coe(2);
out(Ind1) = ( log( (log(mu(Ind1))-Coe(4))/Coe(1) ) - Coe(3)) /Coe(2);
out(Ind2) = log(mu(Ind2));


function out = dlinkExt(mu, Coe)
CC = -1.5;
out = ones(size(mu));
Ind1 = mu <= exp(CC);
Ind2 = mu > exp(CC);

% out(Ind1) = 1 / exp(CC);
out(Ind1) = 1 ./ ( Coe(2) * mu(Ind1) .* (log(mu(Ind1)) - Coe(4)) );
out(Ind2) = 1 ./ mu(Ind2);


function out = ilinkExt(eta, Coe)
CC = -1.5;
out = ones(size(eta));
Ind1 = eta <= CC;
Ind2 = eta > CC;

% out(Ind1) = (eta(Ind1)-CC)*exp(CC) + exp(CC);
out(Ind1) = exp(Coe(1)*exp(Coe(2)*eta(Ind1)+Coe(3)) + Coe(4));
out(Ind2) = exp(eta(Ind2));
























function x = constrain(x,lower,upper)
% Constrain between upper and lower limits, and do not ignore NaN
x(x<lower) = lower;
x(x>upper) = upper;