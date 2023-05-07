% Author: Yu Chen
% Date: Mar 6th 2019 @ CNBC

function [beta, history] = ADMM_glm_lasso(X_train, Y_train, model_list, beta_init, opt)

rho = opt.rho;
zeta = opt.zeta;
lambda = opt.lambda;

for id = 1:length(model_list)
    model_list{id}.lambda = lambda;
end

if lambda == 0
    maxIter = 1;
    beta_loop = 1;
else
    maxIter = opt.maxIter;
    beta_loop = 1;
end
% optimization parameters.
tau_incr = 2;
tau_decr = 0.5;
prim_dual_thr = 10;
eps_abs = 0.001;
eps_rel = 0.001;
eta = 0.999;

% statistical model parameters
B = length(model_list);
d = size(X_train{1}, 2);
nB = size(X_train{1}, 1);
g_list = get_model_par(model_list, 'channel_scalar');

% initialization of variables.
D = build_D(B, d, g_list);
% beta = 10*randn(d, B);
beta = beta_init; 
z = zeros((B-1)*d,1);

% the smart start when lambda = lambda_max 
if opt.smart_start
opt_w.rho = rho;
opt_w.zeta = zeta;
opt_w.B = B;
opt_w.d = d;
w = get_optimal_w(X_train, Y_train, D, beta, opt_w);
else
w = randn((B-1)*d,1);
end

[obj, NLL, penalty] = calc_objective(X_train, Y_train, beta, D, lambda, zeta);
meanLL = -NLL / length(model_list) / length(model_list{1}.Trial_Index_train);
fprintf('%d: obj: %.3f (%.3f, %.3f) \t rho: %.4f \n', ...
    0, obj, meanLL, penalty, rho)

for k = 1:maxIter
zold = z;

%% update beta
opt_Newton.maxiter = 40;
opt_Newton.rho = rho;
opt_Newton.zeta = zeta;
opt_Newton.d = d;
opt_Newton.acc = 1e-4;
opt_Newton.warm_start = 1;

% if k < 5
%     opt_Newton.warm_start = 0;
% else
%     opt_Newton.warm_start = 1;
% end

%+++++ the loop for the beta +++++
if k > 10;  beta_loop = 1;  end

for iii = 1:beta_loop

for i = 1:B
beta(:,i) = ADMM_glm_Newton(X_train{i}, Y_train{i}, D, beta, z, w, opt_Newton, i);
end

%% update z
z = soft_thresholding(D, beta, w, rho, lambda);

end
%% update w
w = w + z - D*beta(:);

%% convergence and residuals analysis

r = z - D*beta(:);
s = rho * D' * (z - zold); 

history.prim_res(k) = norm(r);
history.dual_res(k) = norm(s); 

% keep track of the convergence
[history.obj(k), NLL, penalty] = calc_objective(X_train, Y_train, beta, D, lambda, zeta);

%------- if there is something weird happen, we should restart the calulation
if history.obj(k) > 1e10
    beta = 0*randn(d, B);
    z = zeros((B-1)*d,1);
    w = randn((B-1)*d,1);
    disp('WEIRD thing happens, restart calculation')
end %----------------------------------------

meanLL = -NLL / length(model_list) / length(model_list{1}.Trial_Index_train);

if mod(k, 100) == 0
fprintf('%d: obj: %.3f (%.3f, %.3f) \t pres: %.3f \t dres: %.3f \t rho: %.4f \n', ...
    k, history.obj(end), meanLL, penalty, history.prim_res(end), history.dual_res(end), rho)
end

% relative residuals
eps_prim = sqrt(d)*eps_abs + eps_rel* max( norm(z), norm(D*beta(:)) );
eps_dual = sqrt(B*nB)*eps_abs + eps_rel*rho*norm(w'*D);
history.eps_prim = eps_prim;
history.eps_dual = eps_dual;

% update rho
if k > 25 && opt.adaptive_rho %++++++ update rho after a while %++++++
if history.prim_res(k) < 1e-6 || history.dual_res(k) < 1e-6
    0; % do nothing
elseif history.prim_res(k) > prim_dual_thr * history.dual_res(k)
    rho = rho * tau_incr;
    w = w / tau_incr;
elseif history.dual_res(k) > prim_dual_thr * history.prim_res(k)
    rho = rho * tau_decr;
    w = w / tau_decr;
end
end %++++++ update rho after a while %++++++ 


if history.prim_res(k)<eps_prim && history.dual_res(k)<eps_dual
    break
end

end

%% output
figure('position',[100, 100, 1500, 400])
subplot(131)
plot(history.obj, '-o')
title('objective')
subplot(132)
semilogy(history.prim_res, '-o')
title('primal residual')
subplot(133)
semilogy(history.dual_res, '-o')
title('dual residual')

end

%============================================================================================
function w_star = get_optimal_w(X, Y, D, beta, opt)
B = opt.B;
d = opt.d;
rho = opt.rho;
zeta = opt.zeta;

for i = 1:B
mu = sigmoid(X{i} * beta(:,i));
glm_grad(:,i) = X{i}' * (mu - Y{i}) + zeta*beta(:,i);  
end

w_star = (D*D') \ D * glm_grad(:); 
w_star = w_star / rho; % be very careful about rho * X\b, it is not the same as rho * inv(X)*b
vv = rho * D' * w_star;

end

function b = ADMM_glm_Newton(X, Y, D, beta, z, w, options, i)
d = options.d;
rho = options.rho;
zeta = options.zeta;
acc = options.acc;
Di = D(:,((i-1)*d+1):(d*i));

if options.warm_start
b = beta(:,i);
else
b = zeros(d,1); 
end

obj(1) = -Inf;
%% Newton
Hessian_LagrangianAugment = rho * Di' * Di;

for iter = 1:options.maxiter

mu = sigmoid(X*b);
grad_logistic = X' * (mu - Y) + zeta*b; 
grad_LagrangianAugment = D*beta(:) - z - w; % beta(:) vectorize the beta matrix
grad_LagrangianAugment = rho * Di' * grad_LagrangianAugment;
grad = grad_logistic + grad_LagrangianAugment;

diag = mu .* (1-mu);
Hessian_logistic = X' * (X.*diag) + zeta*eye(d);
Hessian = Hessian_logistic + Hessian_LagrangianAugment;
delta = -Hessian \ grad;
% COND = cond(Hessian, 2); 
% RCOND = rcond(Hessian);

b = b + delta;
beta(:,i) = b;
obj(end+1) = calc_beta_obj(X, Y, beta, D, z, w, i, zeta, rho);

if abs( obj(end) - obj(end-1) ) < acc % 1e-5
    break
end

end %--------- end of the iteration -----------

% fprintf('obj: %d niter: %d  %3.5f\n', i, iter, obj(end))
% figure
% plot(obj)
% title(['Model: ' num2str(i)])

end

function obj = calc_beta_obj(X, Y, beta, D, z, w, i, zeta, rho)
b = beta(:,i);
ll = log_likelihood_logistic(X, Y, b);

L2 = zeta * (b'*b);

largrangin_augmentation = z + w - D*beta(:);
largrangin_augmentation = rho/2 * (largrangin_augmentation' * largrangin_augmentation);

obj = -ll + L2 + largrangin_augmentation;

end

function z = soft_thresholding(D, beta, w, rho, lambda)

x = D * beta(:) - w;

if lambda == 0
    thr = 0;
else
    thr = lambda/rho;
end

z = wthresh(x,'s',thr);
end

function [obj, NLL, penalty] = calc_objective(X, Y, beta, D, lam, zeta)
B = length(X);

NLL = 0;
for i = 1:B
NLL = NLL - log_likelihood_logistic(X{i}, Y{i}, beta(:,i));
end

L2 = zeta * beta(:)' * beta(:);

penalty = lam * norm(D * beta(:), 1);

obj = NLL + L2 + penalty;
end

function D = build_D(B, d, g_list)
g = 1./diff(g_list);

D = sparse((B-1)*d, B*d);

for i = 1:(B-1)
    
    r = ((i-1)*d+1):i*d;
    c = r;
    matind = sub2ind([(B-1)*d, B*d], r, c);
    D(matind) = g(i);

    r = ((i-1)*d+1):i*d;
    c = r + d;
    matind = sub2ind([(B-1)*d, B*d], r, c);
    D(matind) = -g(i);

end

% figure
% imagesc(D)

end

function ll = log_likelihood_logistic(X, y, beta)

eta = X*beta;
second_term = 1 + exp(eta);
second_term = log(second_term);

ll = y' * eta - sum( second_term );

% NLL = -log_likelihood;
% LL = log_likelihood(X,y,beta,'bernoullis', 'offset', offset);

end











