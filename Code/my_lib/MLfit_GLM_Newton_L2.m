function prsML = MLfit_GLM_Newton_L2(X, Y, lam_diag, opts)

nanind = sum(isnan([X Y]), 2) ~= 0;
X(nanind,:) = [];
Y(nanind,:) = [];

lfunc = @(prs)Loss_GLM_logli_exp_Newton_L2(prs, X, Y, lam_diag); % loss function for exponential nonlinearity

prs0 = abs(5*rand(size(X, 2), 1));
% load('E:\CRCNS_Project\prs0.mat')
[prsML, ~] = fminunc(lfunc, prs0, opts); % find ML estimate of params


