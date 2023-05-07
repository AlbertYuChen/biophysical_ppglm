% Author: Yu Chen
% Date: Augus-30-2018 @ MI

function Tau = KS_Haslinger_on_trials(stim, spike_train, KS_test_trials, gg0, KS_valid_point_range)

zeta_list = [];
Tau.U1 = [];
Tau.U2 = [];
Tau.interval_Trial_Start_End = {};
rng(623)
% for ii = 1:size(spike_train_test, 1)
for ii = KS_test_trials
    % Add here for NEURON simulation, where each trial has its own stim. 
    if size(stim,2) == 1
        stim_ind = 1; 
    else
        stim_ind = ii;
    end
    
    p_k_i_list = neuron_forward(gg0, stim(:,stim_ind), spike_train(ii,:));
    q_k_i_list = -log(1-p_k_i_list);
    spk_time = find(spike_train(ii,:));
    valid_spk_time = intersect(spk_time, KS_valid_point_range);
    
    trial_zeta_list = [];
    for spkInd = 2:length(valid_spk_time)
        zeta = q_k_i_list( (valid_spk_time(spkInd-1)+1) : (valid_spk_time(spkInd)-1) );
        zeta = sum(zeta);
        r = rand(1);
        p_last = p_k_i_list(valid_spk_time(spkInd));
        delta = -log(1-r*p_last);
        trial_zeta_list = [trial_zeta_list, zeta+delta];
        % structure of interval_Trial_Start_End
        % trial ID, interval start time, interval end time, interval id, number of intervals
        Tau.interval_Trial_Start_End{end+1} = [ii, (valid_spk_time(spkInd-1)), valid_spk_time(spkInd), ...
                        (spkInd-1), length(valid_spk_time)-1 ]; 
    end
    Tau.U1 = [Tau.U1, trial_zeta_list(1:end-1)];
    Tau.U2 = [Tau.U2, trial_zeta_list(2:end)];
    zeta_list = [zeta_list, trial_zeta_list];
end


%%
Tau.zeta_list = zeta_list;
Tau.Uniform_list = 1- exp(-zeta_list);

[eCDF, zvals] = ecdf(zeta_list);
mCDF = 1-exp(-zvals); %Model CDF at z values.

figure('Position', [200, 200, 1400, 300])
subplot(141)
plot(mCDF, eCDF, 'Linewidth', 3)
hold on	 

% 0.95% CI
plot([0 1], [0 1]+1.36/sqrt(length(zeta_list)), '--', 'Color', [.7 .7 .7])		%Upper confidence bound.
plot([0 1], [0 1]-1.36/sqrt(length(zeta_list)), '--', 'Color', [.7 .7 .7])		%Lower confidence bound.

% 0.999% CI
% plot([0 1], [0 1]+1.95/sqrt(length(Tau.Zscr)), '--', 'Color', [.7 .7 .7])		%Upper confidence bound.
% plot([0 1], [0 1]-1.95/sqrt(length(Tau.Zscr)), '--', 'Color', [.7 .7 .7])		%Lower confidence bound.

xlabel('Model CDF')
ylabel('Empirical CDF')
axis([0 1 0 1]);
title('KS test')
grid on; set(gca,'FontSize', 12)

% KS test p-value
epsilon95 = 1.36/sqrt(length(zeta_list));
KS_p_value95 = 2*exp(-2*length(zeta_list)*epsilon95^2);
epsilon = max(abs(mCDF - eCDF));
KS_p_value = 2*exp(-2*length(zeta_list)*epsilon^2);

disp(['epsilon  95%: ' num2str(epsilon95) '  p-value: ' num2str(KS_p_value95)])
disp(['epsilon data: ' num2str(epsilon) '  p-value: ' num2str(KS_p_value)])

%--- iid test ---

Ui1 = 1 - exp(-Tau.U1);
Ui2 = 1 - exp(-Tau.U2);
subplot(142)
plot(Ui1, Ui2, '.', 'Markersize', 4)
xlabel('U_i')
ylabel('U_{i+1}')
title('IID test')
% axis equal; 
axis([0 1 0 1])
set(gca,'FontSize', 12)

%------ histogram -------
subplot(143)
t = 0:0.01:10;
exppdf = exp(-t);
plot(t, exppdf, 'Color', [0 0.7 0], 'Linewidth', 2)
hold on
histogram(zeta_list, 'Normalization', 'pdf', 'FaceColor', [0 0.4470 0.7410], 'BinWidth', 0.1)
xlim([0 5])
legend('Ideal curve')
title('PDF of \tau')
grid on; set(gca, 'FontSize',12)

subplot(144)
plot([0 1],[1 1], 'Color', [0 0.7 0], 'Linewidth', 2)
hold on
histogram(1- exp(-zeta_list), 'Normalization', 'pdf', 'FaceColor', [0 0.4470 0.7410], 'BinWidth', 0.02)
legend('Ideal curve', 'location', 'southeast')
title('PDF of 1 - exp\{-\tau\}')
xlim([0 1])
grid on; set(gca, 'FontSize',12)

end

function output_p = neuron_forward(gg0, stim, single_spike_train)
output_p = zeros(1, length(stim));
if size(gg0.k,2)==1 && size(gg0.h,2)==1
    output_p = neuron_simple_forward(gg0, stim, single_spike_train);
else
    output_p(gg0.window_timepoint) = neuron_time_varying_forward(gg0, stim, single_spike_train);
end
end

function output_p = neuron_time_varying_forward(gg0, stim, single_spike_train)
k_response = stim_conv(stim, gg0.k, gg0.k_delay);
h_response = hist_conv(single_spike_train', gg0.h);

eta = k_response + h_response + gg0.dc;
output_p = sigmoid(eta, false);

linear_index = sub2ind(size(output_p), gg0.window_timepoint, gg0.window_timepoint);
output_p = output_p(linear_index);
end


function output_p = neuron_simple_forward(gg0, stim, single_spike_train)

if isfield(gg0, 'k_delay')
    k_delay = gg0.k_delay;
else
    k_delay = 0;
end

k_response = stim_conv(stim, gg0.k, k_delay);
h_response = hist_conv(single_spike_train', gg0.h);

eta = k_response + h_response + gg0.dc;
output_p = sigmoid(eta, false);
end




































