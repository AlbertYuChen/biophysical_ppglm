function [logli, dL, H] = Loss_GLM_logli_exp_Newton(prs, X, Y)
dt = 0.001;
Itot = X * prs; 
rr = exp(Itot);
% ---------  Compute log-likelihood ---------------------------------
% logli = sum(rr)*dt - sum(Itot(Y));
logli = sum(rr)*dt - Itot'*Y;
% Combine terms
dL = (rr'*X)'*dt - (Y'*X)';

% ---------  Compute Hessian -----------------
rlen = length(rr);
ddrrdiag = spdiags(rr,0,rlen,rlen); 
H = (X'*ddrrdiag*X)*dt; 













