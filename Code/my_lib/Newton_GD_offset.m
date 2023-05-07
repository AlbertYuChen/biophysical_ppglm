% Author: Yu Chen
% Date: Jul 25th 2018 @ Beijing

% X = X_const;
X = Xbl;
C = baseline_const;
y = gg0.sps(irt_ind_cum);
beta = zeros(size(X,2), 1);

N = length(X);
%% Gradient descent
% lr = 0.5;
% 
% for i = 1:1
%     
% sig = sigmoid(X*beta);
% grad = -X' * (y - sig) / 2001; 
% beta = beta - lr*grad;
% end

%% Newton
for i = 1:20

sig = sigmoid(X*beta+C);
grad = -X' * (y - sig); 

diag = sig .* (1-sig);
Hessian = X' * (X.*diag);

delta = Hessian \ grad;
COND = cond(Hessian, 2); 
RCOND = rcond(Hessian);
fprintf('Iter: %d \t COND: %3e \n', i, COND);

beta = beta - delta;
end

%%
% figure 
% plot(prsML)
% hold on
% plot(beta)
% 
% figure
% plot(prsML - beta)
%% log-likelihood
% beta = ones(21,1)*0.1;
% eta = X*prsML0;
eta = X*beta+C;
second_term = 1 + exp(eta);
second_term = log(second_term);

log_likelihood = y' * eta - sum( second_term );
NLL = -log_likelihood / N;

%%
gg0.blw = beta;
gg0.bl = gg0.BLBbas*gg0.blw;

figure('Position', [300, 300, 400, 250]);
plot(gg0.blt, gg0.bl, 'r');
% ylim([-15 5])
% xlim([0 gg0.iht(end)])
grid on


















