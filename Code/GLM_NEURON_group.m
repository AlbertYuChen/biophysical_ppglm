% Author: Yu Chen
% Date: Feb 23rd 2019 @ CNBC

close all;
% clearvars -except seed_id dateOut   dateOut_list    global_id   project_workspace_path ...
%                   channelName   dataSetName    channelName_list   neuron_id_list   N_list; 
%%
%-------------- bhalla --------------
datasetName = 'bhalla_scaled_9520_downsampled'; 
channelName_list = {'kA', 'kca3', 'kfasttab', 'kslowtab', 'lcafixed', 'nafast'};
dateOut_list = {'_Mar11_', '_Mar26_', '_Mar26_', '_Mar20_', '_Mar23_', '_Mar27_'};

N_list{1} = [1 2 3 4 5 8 6 7 9 10];          % kA
N_list{2} = [11 12 13 14 15 18 16 17 19 20]; % kca3
N_list{3} = [25 28 26 27 29 30];             % kfasttab
N_list{4} = [31 32 33 34 35 38 36 37 39 40]; % kslowtab
N_list{5} = [41 42 43 44 45 48 46 47 49 50]; % lcafixed
N_list{6} = [54 55 58 56 57 59 60];          % nafast

%-------------- alon --------------
% datasetName = 'alon_scaled_8015_downsampled'; 
% channelName_list = {'bk', 'cah', 'car', 'iA', 'iH', 'kslow', 'na', 'sk'};
% dateOut_list = {'_Mar28_', '_Mar28_', '_Mar28_', '_Mar28_', '_Mar28_', '_Mar28_', '_Mar28_', '_Mar28_'};
% 
% N_list{1} = [1 2 3 4 5 8 6 7 9];             % bk
% N_list{2} = [11 12 13 14 15 18 16 17];       % cah
% N_list{3} = [21 22 23 24 25 28 26 27];       % car
% N_list{4} = [31 32 33 34 35 38];             % iA
% N_list{5} = [43 44 45 48 46 47 49];          % iH
% N_list{6} = [51 52 53 54 55 58 56 57];       % kslow
% N_list{7} = [66 67 68];                      % na
% N_list{8} = [71 72 73 74 75 78 76 77 79];    % sk

%-------------- simulation --------------
% datasetName = 'channel_scaled_9520_simulation'; 
% datasetName = 'rand_channel_scaled_9520_simulation_11'; 
% seed_id = 12;
% disp(['Training seed_id ', num2str(seed_id)])
% dataSetName = ['rand_channel_scaled_9520_simulation_', num2str(seed_id)];
% channelName_list = {'sim1'};
% dateOut_list = {'_Apr2_'};
% N_list{1} = 1:10;             % sim1 based on ballon kca3

% -------------------------
channel_id = 1;
channelName = channelName_list{channel_id};
dateOut = dateOut_list{channel_id};
neuron_id_list = N_list{channel_id};
CV_fold = 1;
global_id = 1;
spikes_data_file = [project_workspace_path 'Data/simulation_study/'...
    dataSetName '/' dataSetName '.mat'];
outputFolder = [project_workspace_path, 'Output/simulation_study/'...
    dataSetName '_GLM/' channelName dateOut num2str(global_id) '/']; 
mkdir(outputFolder)
%% check working environment
% I'm swtching between PC, Mac, and CNBC cluster.
project_workspace_path = Initialization_Env;
%% load dataset

dt = 0.001;
valid_point_range = (601:3000)';
model_list = {};
for neuron_index = neuron_id_list
model.neuron_index = neuron_index;
model.dt = dt; 
model.filepath = spikes_data_file;
model.valid_point_range = valid_point_range;
model_list{end+1} = model;
end
[all_spike_train, stim_list, model_list, ISIs] = load_spike_train(model_list);

% this is redundant, but easy to check MC 
load(spikes_data_file); 

%% select train, test sets
rng(global_id); % fix the shuffle for the perputation. 

% Trial_Index_all_rand = model.Trial_Index_all;
% N_train_trials = round(length(model.Trial_Index_all) * 2/3 );

% Trial_Index_all_rand = model_list{1}.Trial_Index_all;
Trial_Index_all = model_list{1}.Trial_Index_all;
Trial_Index_all_rand = Trial_Index_all( randperm(length(Trial_Index_all)) );

% N_train_trials = round(100 * 2/3 );
N_train_trials = 67;
Trial_Index_train = Trial_Index_all_rand(1:N_train_trials); 
Trial_Index_test = Trial_Index_all_rand((N_train_trials+1):end); 

%% preview data 
% probe = 10;
% visualization_spiketrains(model_list{probe}.Trial_Index_all, all_spike_train{probe}, ISIs{probe}, [0 3.2])
% 
% figure
% stem(sum(all_spike_train{probe},2));
% ylim([50, 200])

%% X, Y build
for id = 1:length(model_list)
%-------------- Train dataset. --------------
%+++ Xstim_list %+++
% alon_scaled_8015_downsampled_GLM  50 + 10
% bhalla_scaled_9520_downsampled 150 + 10
% channel_scaled_9520_simulation 150 + 10
model_list{id}.nkt = 150;  
model_list{id}.nk = 10;

model_list{id}.Trial_Index_train = Trial_Index_train;

[Xstim, model_list{id}] = build_X_stim(stim_list{id}, Trial_Index_train, model_list{id});
%+++ Xsp %+++
% alon_scaled_8015_downsampled_GLM  130 + 10
% bhalla_scaled_9520_downsampled 150 + 10
% channel_scaled_9520_simulation 150 + 10
model_list{id}.nht = 150; 
model_list{id}.nh = 10;

[Xsp, model_list{id}] = build_X_history(all_spike_train{id}, Trial_Index_train, model_list{id});
% Xdc 
model_list{id}.ndc = 1;

Xdc = ones(length(model_list{id}.valid_point_range)*length(Trial_Index_train), model_list{id}.ndc);

model_list{id}.dim = model_list{id}.nk + model_list{id}.ndc + model_list{id}.nh;
% X, Y
Y = all_spike_train{id}(Trial_Index_train, model_list{id}.valid_point_range);

Y_list{id} = reshape(Y', [], 1);
X_list{id} = [Xstim, Xdc, Xsp];


%-------------- Validation dataset. --------------
% Xstim_list
Trial_Index_test = setdiff(model_list{id}.Trial_Index_all, model_list{id}.Trial_Index_train);
[Xstim, model_list{id}] = build_X_stim(stim_list{id}, Trial_Index_test, model_list{id});
% Xsp
[Xsp, model_list{id}] = build_X_history(all_spike_train{id}, Trial_Index_test, model_list{id});
% Xdc 
Xdc = ones(length(model_list{id}.valid_point_range)*length(Trial_Index_test), model_list{id}.ndc);

model_list{id}.dim = model_list{id}.nk + model_list{id}.ndc + model_list{id}.nh;
% X, Y
Y = all_spike_train{id}(Trial_Index_test, model_list{id}.valid_point_range);

Y_list_test{id} = reshape(Y', [], 1);
X_list_test{id} = [Xstim, Xdc, Xsp];

end

%% joint uniform X, Y estimation
% The data is from getLambdaList.m
model_0_file = [outputFolder, model_list{1}.channel_type, '_0', '.mat'];
if exist( model_0_file, 'file')
disp('load uniform model...')
load([outputFolder, model_list{1}.channel_type, '_0', '.mat']);
betaUni = beta;
else
disp('0 model does not exist')
end
%% get lambda_max
% The data is from getLambdaList.m
if exist([outputFolder, 'lambda_list.mat'], 'file')
    load([outputFolder, 'lambda_list.mat'])
else
    lamMax = get_lambdaMax(X_list, Y_list, betaUni, model_list);
    lamMax = lamMax * 2;
    lambda_list = [0, exp((log(lamMax)-20):1:log(lamMax))];
    save([outputFolder, 'lambda_list.mat'], 'lambda_list')
end

%% start ADMM
% for lambda_id = 1:length(lamda_list)
beta = 0*randn(model_list{1}.dim, length(model_list));

for lambda_id = 1:22

lambda = lambda_list(lambda_id);
opt.lambda = lambda;
opt.zeta = 0.001;
opt.maxIter = 12345;

if lambda_id <= 19
opt.adaptive_rho = 1;
opt.rho = opt.lambda/10;
opt.smart_start = 0;
else % for last few trials 20 21 22
beta = betaUni;
opt.adaptive_rho = 1;
opt.rho = opt.lambda/15;
opt.smart_start = 1;
end

B = length(model_list);
d = model_list{1}.dim;

tStart = tic;

[beta, history] = ADMM_glm_lasso(X_list, Y_list, model_list, beta, opt);

tEnd = toc(tStart);

fprintf('ellapse time: %d min %d s\n', floor(tEnd/60), floor(mod(tEnd,60)) );

% assemble results into model_list, and plot results.
model_list = assemble_filters(beta, model_list);

% A post analysis is in script GLM_CrossValidation.m
for id = 1:length(model_list)
    model_list{id}.log_likelihood_train = log_likelihood(...
        X_list{id}, Y_list{id}, beta(:,id), 'bernoulli');
    model_list{id}.log_likelihood_test = log_likelihood(...
        X_list_test{id}, Y_list_test{id}, beta(:,id), 'bernoulli');
end

%% Save models.
% plot_models(model_list)
% plot_training_history(history)

fprintf('save checkpoint: %d \n', lambda_id)
save([outputFolder, model_list{1}.channel_type, '_', num2str(lambda_id), '.mat'], ...
        'model_list', 'beta', 'lambda', 'history')

end















