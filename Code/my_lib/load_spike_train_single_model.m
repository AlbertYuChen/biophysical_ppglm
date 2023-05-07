% Author: Yu Chen
% Date: Feb 23rd 2019 @ CNBC

function [all_spike_train, stim, model, ISIs] = load_spike_train_single_model(model)
neuron_index = model.neuron_index;
load(model.filepath);

if exist('MC2', 'var')
    MC = MC2;
else
    disp('MC2 does not exist')
end

if ~exist('stim', 'var')
    disp('Stimulus variable does not exist')
end


if isfield(MC, 'dataset')
% This is for Nate's dataset name e.g.  spiketimes_kA_1.5.csv
experimentName = MC(neuron_index).dataset;
if strcmp(experimentName((end-3):end), '.csv')
experimentName = experimentName(1:(end-4));
end
C = strsplit(experimentName, '_');
channel_type = C{2};
channel_scalar = str2double( C{3} );

model.experimentName = experimentName;
model.channel_type = channel_type;
model.channel_scalar = channel_scalar;
end

fprintf('Load neuron: %d \t %s \n', neuron_index, experimentName);

% ------------------

dt = model.dt;
stim_len = length(stim);
stim_time = dt:dt:stim_len*dt;

model.stim_len = stim_len;
model.stim_time = stim_time;

% figure
% plot(stim_time, stim)

if isfield(model, 'valid_point_range')
    valid_point_range = model.valid_point_range;
else
    valid_point_range = (1:stim_len)';
end

Trial_Index_all = 1:max(MC(neuron_index).spikeIndices);
model.Trial_Index_all = Trial_Index_all;

%% preprocessing
all_spike_train = []; 
ISIs = [];

for ii = 1:max(MC(neuron_index).spikeIndices)
    spike_index_array = MC(neuron_index).spikeIndices == ii;
    single_trial_spike_time = MC(neuron_index).spikeTimes(spike_index_array);
    single_trial_spike_time = round(single_trial_spike_time);
    single_trial_spike_time(single_trial_spike_time==0) = []; % this is first found in alon 'spiketimes_na_1.5.csv' 
    single_trial_spike_time(single_trial_spike_time>valid_point_range(end)) = []; % this is for Tripthy data.
    
    x = zeros(1, stim_len);
    x(single_trial_spike_time) = 1;
    all_spike_train = [all_spike_train; x]; 
    
    single_trial_spike_time_dt = single_trial_spike_time * dt;
    ISIs = [ISIs diff(single_trial_spike_time_dt)]; 
end

end
