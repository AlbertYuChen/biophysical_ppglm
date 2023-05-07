% Author: Yu Chen
% Date: June 2nd 2018 @ CNBC CMU

function visualization_spiketrains(Trial_Index, all_spike_train, ISIs, customized_range)
spike_rec_res = 0.001;
stim_rec_res = 0.001;
stim_len = size(all_spike_train,2);
stim_time = 0.001:stim_rec_res:stim_len*stim_rec_res;


%%-------------- PSTH -----------------
PSHT_time_bin_width = 0.001;
PSHT_time_bin = PSHT_time_bin_width:PSHT_time_bin_width:stim_len*stim_rec_res;
x_lim_1 = [0 1];
x_lim_2 = [1 2];

PSTH_reshape1 = plot_PSTH(all_spike_train, PSHT_time_bin, stim_time, [], x_lim_1, 2);
% figure
% plot(PSHT_time_bin, PSTH_reshape1 ,'Color',[0 0.5 0], 'LineWidth', 1)
% title('PSTH')


%-------------- Raster plot -----------------
% figure('Position', [300, 300, 800, 500]);
% imagesc(stim_time, Trial_Index, all_spike_train) 
% colormap(1-gray); 

%% ------------ Raw data raster ------------
figure('Position', [100, 100, 1800, 800]);
subplot(211)
% imagesc(stim_time, Trial_Index, all_spike_train)			
% colormap(1-gray);	
for ii = 1:size(all_spike_train,1)
plot(find(all_spike_train(ii,:))*spike_rec_res, Trial_Index(ii)*ones(sum(all_spike_train(ii,:)),1),...
    's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'Markersize', 2); hold on
end
xlim(x_lim_1)
grid on; set(gca,'FontSize',12)
%------------------------------------
subplot(212)
plot(PSHT_time_bin, PSTH_reshape1 ,'Color',[0 0.5 0], 'LineWidth', 1)

% yyaxis right
% plot(stim_time, stim)

xlim(x_lim_1)

%%
figure('Position', [100, 100, 1800, 800]);
subplot(211)
for ii = 1:size(all_spike_train,1)
plot(find(all_spike_train(ii,:))*spike_rec_res, Trial_Index(ii)*ones(sum(all_spike_train(ii,:)),1),...
    's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'Markersize', 2); hold on
end
xlim(x_lim_2)
grid on; set(gca,'FontSize',12)
%------------------------------------
subplot(212)
plot(PSHT_time_bin, PSTH_reshape1 ,'Color',[0 0.5 0], 'LineWidth', 1)
% yyaxis right
% plot(stim_time, stim)
xlim(x_lim_2)

%%
if nargin >= 3
% ---------- ISI ------------
% ISI_bin_width = 0.4; % Human
% ISI_bin_width = 0.005;
figure
histogram(ISIs, 40, 'Normalization', 'pdf');	%...compute histogram of ISIs,
end


%%
if nargin >= 4
figure('Position', [100, 100, 1800, 800]);
subplot(211)
for ii = 1:size(all_spike_train,1)
plot(find(all_spike_train(ii,:))*spike_rec_res, Trial_Index(ii)*ones(sum(all_spike_train(ii,:)),1),...
    's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'Markersize', 2); hold on
end
xlim(customized_range)
ylabel('Trial index')
xlabel('Time [s]')
grid on; set(gca,'FontSize',16)
%------------------------------------
subplot(212)
plot(PSHT_time_bin, PSTH_reshape1 ,'Color',[0 0.5 0], 'LineWidth', 1)
% yyaxis right
% plot(stim_time, stim)
xlim(customized_range)
ylabel('Firing rate [spk/ms]')
xlabel('Time [s]')
grid on; set(gca,'FontSize',16)
end




















