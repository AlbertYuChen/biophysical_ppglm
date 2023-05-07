% Author: Yu Chen
% Date: June 3rd 2018 @ CNBC
% The algorithm is decribed in Citi et al. 2014 Appendix C

function beta = Citi_Newton(X, Y)
max_iteration = 50;

[n, p] = size(X);
Yrho = 1 - Y/2; 

tiny = realmin('double')^.25; % keep fourth powers from under/overflowing
lowerBnd = log(tiny); 
upperBnd = -lowerBnd;
ilink = @(eta) exp(clip(eta,lowerBnd,upperBnd));

%% initialize beta
beta = zeros(p, 1);

%% begin iterations

for iter = 1:max_iteration
rho_lambda = Yrho .* ilink( X * beta );
Gradient = X' * (Y - rho_lambda);
Hessian = -X' * (rho_lambda .* X);
beta = beta - Hessian \ Gradient;
end


% link = @(mu) log(mu);
% dlink = @(mu) 1 ./ mu;
% % keep mu = ilink(eta) in [tiny, 1/tiny];

%-----------------------------------------------------------------------
function x = clip(x,lower,upper)
% Constrain between upper and lower limits, and do not ignore NaN
x(x<lower) = lower;
x(x>upper) = upper;












