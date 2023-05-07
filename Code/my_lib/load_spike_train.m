% Author: Yu Chen
% Date: Feb 23rd 2019 @ CNBC


function [all_spike_train, stim, model_new, ISIs] = load_spike_train(model)

all_spike_train = {};
stim = {};
model_new = {};
ISIs = {};

for i = 1:length(model)
    [all_spike_train{i}, stim{i}, model_new{i}, ISIs{i}] = load_spike_train_single_model(model{i});
end

end






