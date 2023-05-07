% Author: Yu Chen
% Date: Aug 27th 2017 @ CNBC

function Itot = Itot_function_IMI(glmprs, Stim, spktr)

slen = size(Stim,1);

lastspk = 0;
last2spk = 0;
last3spk = 0;
last4spk = 0;

Istm = glmprs.kdc;
Itot = Istm; 

%% --------------- Append zeros to filters ---------------
if isfield(glmprs, 'ih1')
    glmprs.ih1 = [glmprs.ih1; zeros(300, 1)];
end
if isfield(glmprs, 'ih2')
    glmprs.ih2 = [glmprs.ih2; zeros(300, 1)];
end
if isfield(glmprs, 'ih3')
    glmprs.ih3 = [glmprs.ih3; zeros(300, 1)];
end
if isfield(glmprs, 'ih4')
    glmprs.ih4 = [glmprs.ih4; zeros(300, 1)];
end
if isfield(glmprs, 'ih1h1h1')
    glmprs.ih1h1h1 = [glmprs.ih1h1h1; zeros(800, 1)];
end

% --------------- Run dynamics ---------------
for ispk = [find(spktr) slen]
    
    %---- update synthesised firing rate ----
    iiPostSpk = (lastspk+1):ispk;        
    if isfield(glmprs, 'ih3')
        Itot(iiPostSpk) = Itot(iiPostSpk) + glmprs.ih1(iiPostSpk - lastspk)... 
            + glmprs.ih2(iiPostSpk - last2spk) + glmprs.ih3(iiPostSpk - last3spk);
    else
        disp('does not contain history dpendency function')
    end
    % --  Update # of samples per iter ---
    last4spk = last3spk;
    last3spk = last2spk;
    last2spk = lastspk;
    lastspk = ispk;

end















