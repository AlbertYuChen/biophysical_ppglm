% Author: Yu Chen
% Date: Aug 27th 2017 @ CNBC

function ItotNLN = Itot_function_NLNlink(glmprs, Stim, spktr)

slen = size(Stim,1);
hlen = size(glmprs(1).ih,1); % length of post-spike filter
Istm = glmprs.kdc;
Itot = Istm; 

% --------------- Run dynamics ---------------
for ispk = find(spktr)

    %---- update synthesised firing rate ----
    iiPostSpk = (ispk+1):min(slen, ispk+hlen);     
    Itot(iiPostSpk) = Itot(iiPostSpk)+glmprs.ih(iiPostSpk-ispk);

end


%% NOTE:
% the nonlinear function will be handled in @Tau_cross_validation_NLNlink
% so we just skip the nonlinear steps here. 

ItotNLN = Itot;
% validIndex = Itot<=glmprs.BasisEnd & Itot>=glmprs.BasisBegin;
ItotNLN(ItotNLN>glmprs.BasisEnd) = 0;
ItotNLN(ItotNLN<glmprs.BasisBegin) = -30;    
ItotBasisIndex = round( (Itot - glmprs.BasisBegin) / glmprs.dx);    
ItotNLN = glmprs.f_Link_NLN(ItotBasisIndex);















