% Author: Yu Chen
% Date: June 2nd 2018 @ CNBC CMU

function visualization_spiketrains_DEBUG(Trial_Index, all_spike_train, Tau, outlier_Index)
spike_rec_res = 0.001;
stim_rec_res = 0.001;
stim_len = size(all_spike_train,2);
stim_time = 0.001:stim_rec_res:stim_len*stim_rec_res;


%%-------------- PSTH -----------------
PSHT_time_bin_width = 0.001;
PSHT_time_bin = PSHT_time_bin_width:PSHT_time_bin_width:stim_len*stim_rec_res;
x_lim_1 = [0 1];
x_lim_2 = [0.9 2];

PSTH_reshape1 = plot_PSTH(all_spike_train, PSHT_time_bin, stim_time, [], x_lim_1, 0);
% figure
% plot(PSHT_time_bin, PSTH_reshape1 ,'Color',[0 0.5 0], 'LineWidth', 1)
% title('PSTH')


%-------------- Raster plot -----------------
% figure('Position', [300, 300, 800, 500]);
% imagesc(stim_time, Trial_Index, all_spike_train) 
% colormap(1-gray); 


%% ------------ Raw data raster ------------
% figure('Position', [100, 100, 1800, 800]);
% subplot(211)
% % imagesc(stim_time, Trial_Index, all_spike_train)			
% % colormap(1-gray);	
% for ii = 1:size(all_spike_train,1)
% plot(find(all_spike_train(ii,:))*spike_rec_res, Trial_Index(ii)*ones(sum(all_spike_train(ii,:)),1),...
%     's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'Markersize', 2); hold on
% end
% 
% for jj = outlier_Index
% interval_Trial_Start_End = Tau.interval_Trial_Start_End{jj};
% plot([interval_Trial_Start_End(2)*spike_rec_res; interval_Trial_Start_End(3)*spike_rec_res],...
%     [interval_Trial_Start_End(1); interval_Trial_Start_End(1)], 'r')
% hold on
% end
% 
% xlim(x_lim_1)
% ylabel('Trial index') 
% grid on; set(gca,'FontSize',18)
% %------------------------------------
% subplot(212)
% plot(PSHT_time_bin, PSTH_reshape1 ,'Color',[0 0.5 0], 'LineWidth', 1)
% 
% % yyaxis right
% % plot(stim_time, stim)
% 
% xlim(x_lim_1)
% ylabel('PSTH firing rate [Hz]')
% xlabel('Time [sec]')
% grid on; set(gca,'FontSize',18)
% %%
% figure('Position', [100, 100, 1800, 800]);
% subplot(211)
% for ii = 1:size(all_spike_train,1)
% plot(find(all_spike_train(ii,:))*spike_rec_res, Trial_Index(ii)*ones(sum(all_spike_train(ii,:)),1),...
%     's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'Markersize', 2); hold on
% end
% 
% for jj = outlier_Index
% interval_Trial_Start_End = Tau.interval_Trial_Start_End{jj};
% plot([interval_Trial_Start_End(2)*spike_rec_res, interval_Trial_Start_End(3)*spike_rec_res],...
%     [interval_Trial_Start_End(1) interval_Trial_Start_End(1)], 'r')
% hold on
% end
% 
% 
% xlim(x_lim_2)
% ylabel('Trial index') 
% grid on; set(gca,'FontSize',18)
% %------------------------------------
% subplot(212)
% plot(PSHT_time_bin, PSTH_reshape1 ,'Color',[0 0.5 0], 'LineWidth', 1)
% % yyaxis right
% % plot(stim_time, stim)
% xlim(x_lim_2)
% ylabel('PSTH firing rate [Hz]')
% xlabel('Time [sec]')
% grid on; set(gca,'FontSize',18)

%% ------------ Raw data raster ------------
figure('Position', [100, 100, 1800, 800]);
subplot(211)
% imagesc(stim_time, Trial_Index, all_spike_train)			
% colormap(1-gray);	
for ii = 1:size(all_spike_train,1)
plot(find(all_spike_train(ii,:))*spike_rec_res, Trial_Index(ii)*ones(sum(all_spike_train(ii,:)),1),...
    's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'Markersize', 2); hold on
end

for jj = outlier_Index
interval_Trial_Start_End = Tau.interval_Trial_Start_End{jj};
plot([interval_Trial_Start_End(2)*spike_rec_res; interval_Trial_Start_End(3)*spike_rec_res],...
    [interval_Trial_Start_End(1); interval_Trial_Start_End(1)], 'r')
hold on
end

xlim(x_lim_1)
ylabel('Trial index') 
grid on; set(gca,'FontSize',18)
%------------------------------------
subplot(212)

for ii = 1:size(all_spike_train,1)
plot(find(all_spike_train(ii,:))*spike_rec_res, Trial_Index(ii)*ones(sum(all_spike_train(ii,:)),1),...
    's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'Markersize', 2); hold on
end

for jj = outlier_Index
interval_Trial_Start_End = Tau.interval_Trial_Start_End{jj};
plot([interval_Trial_Start_End(2)*spike_rec_res, interval_Trial_Start_End(3)*spike_rec_res],...
    [interval_Trial_Start_End(1) interval_Trial_Start_End(1)], 'r')
hold on
end


xlim(x_lim_2)
ylabel('Trial index') 
xlabel('Time [sec]')
grid on; set(gca,'FontSize',18) 
















