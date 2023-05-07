function prsML = MLfit_GLM_Newton(X, Y)

lfunc = @(prs)Loss_GLM_logli_exp_Newton(prs, X, Y); % loss function for exponential nonlinearity

[~, Xwidth] = size(X);
prs0 = zeros(Xwidth, 1);
opts = optimset('Gradobj','on','Hessian','on','display', 'iter', 'maxiter', 1000);
% load('E:\CRCNS_Project\prs0.mat')
[prsML, ~] = fminunc(lfunc, prs0, opts); % find ML estimate of params
