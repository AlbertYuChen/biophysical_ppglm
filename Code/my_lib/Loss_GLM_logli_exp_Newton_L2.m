function [logli, dL, H] = Loss_GLM_logli_exp_Newton_L2(prs, X, Y, lam_diag)
% dt = 0.001;
dt = 1;
Itot = X * prs; 
rr = exp(Itot);
% ---------  Compute log-likelihood ---------------------------------
logli = sum(rr)*dt - Itot'*Y + prs'*lam_diag*prs;
% Combine terms
dL = (rr'*X)'*dt - (Y'*X)' + 2*dt*lam_diag*prs;
% ---------  Compute Hessian -----------------
rlen = length(rr);
ddrrdiag = spdiags(rr,0,rlen,rlen); 
H = (X'*ddrrdiag*X)*dt + 2*lam_diag*dt; 













