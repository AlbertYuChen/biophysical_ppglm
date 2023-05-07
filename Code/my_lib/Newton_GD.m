% Author: Yu Chen
% Date: Jul 25th 2018 @ Beijing

X = X_design;
y = gg0.sps(irt_ind_cum);
beta = ones(size(X_design,2),1) * 0;
% beta = [1.8020, -14.0592,   6.3536,  -0.1072,   0.1437,  -0.4509, ...
%            0.5527,  -0.4453,   0.2188,  -0.1460,  -2.3305, -13.6198, ...
%          -19.8999,  -6.9952, -31.5315, -13.3863,  -4.0088,  -4.1005, ...
%           -0.8276,  -2.2776,   1.5278]';

N = length(X);
%% Newton
for i = 1:20

sig = sigmoid(X*beta);
grad = -X' * (y - sig); 

diag = sig .* (1-sig);
Hessian = X' * (X.*diag);

delta = Hessian \ grad;
COND = cond(Hessian, 2); 
RCOND = rcond(Hessian);
fprintf('Condition Number: %3e \n', COND);

beta = beta - delta;
end

%% log-likelihood
eta = X*beta;
second_term = 1 + exp(eta);
second_term = log(second_term);

log_likelihood = y' * eta - sum( second_term );
NLL = -log_likelihood / N;

%%
kw = beta(1:10);
k = gg0.ktbas*kw;

figure('Position', [300, 300, 400, 250]);
plot(k, 'r'); 

%---------------------------- H ---------------------------
hw = beta(12:end);
h = gg0.ihbas*hw; 

figure('Position', [300, 300, 400, 250]);
plot(h,'r'); 






















