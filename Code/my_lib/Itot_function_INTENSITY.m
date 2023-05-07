% Author: Yu Chen
% Date: Aug 27th 2017 @ CNBC

function Itot = Itot_function_INTENSITY(glmprs, Stim, spktr)

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















