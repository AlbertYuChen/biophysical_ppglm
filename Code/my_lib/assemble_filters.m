% Author: Yu Chen
% Date: Feb 24th 2019 @ CNBC

function model_list = assemble_filters(beta_mat, model_list)

if length(model_list) == 1; model_list = {model_list}; end

for i = 1:length(model_list)
    
model = model_list{i};
par = beta_mat(:,i);

nk = model.nk;
nh = model.nh;
ndc = model.ndc;

model.kw = par(1:nk);
model.k = model.kbas*model.kw;

if ndc == 0
    model.dc = 0;
else
    model.dc = par(nk+1);
end

% model.kdc = stim_conv(stim_list{i}, model.k, model.k_delay) + model.dc;

model.hw = par(nk+ndc+1:nk+nh+ndc);
model.h = model.hbas*model.hw;

model_list{i} = model;

end

if length(model_list) == 1; model_list = model_list{1}; end


