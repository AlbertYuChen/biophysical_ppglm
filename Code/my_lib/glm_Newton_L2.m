% Author: Yu Chen
% Date: Aug 11th 2018 @ Beijing

function [beta, ll] = glm_Newton_L2(X,y,opt)

zeta = opt.zeta;
d = size(X,2);

beta = ones(d,1) * 0;
N = length(y);

%% Newton
obj(1) = Inf;

for i = 1:opt.maxiter
mu = sigmoid(X*beta);
grad = -X' * (y - mu) + zeta*beta; 
diag = mu .* (1-mu);
Hessian = X' * (X.*diag) + zeta*eye(d);

delta = Hessian \ grad;

if opt.verbose
[~, MSGID] = lastwarn();
warning('off', MSGID)
end
% COND = cond(Hessian, 2); 
% RCOND = rcond(Hessian);

beta = beta - delta;

obj(end+1) = calc_obj(X, y, beta, zeta);
if abs( obj(end) - obj(end-1) ) < 1e-9
    break
end

end

ll = obj(end);

% figure
% plot(obj, '-o')

end

%% log-likelihood for logistic (other types of distribution have not been implemented yet. )

function obj = calc_obj(X, y, beta, zeta)
N = length(y);
eta = X*beta;
second_term = 1 + exp(eta);
second_term = log(second_term);

log_likelihood = y' * eta - sum( second_term );
obj = log_likelihood + zeta*beta'*beta;
obj = obj / N;
% LL = log_likelihood(X,y,beta,'bernoullis', 'offset', offset);

end




















