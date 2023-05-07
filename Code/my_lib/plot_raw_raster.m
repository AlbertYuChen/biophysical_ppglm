% Author: Yu Chen
% Date: April 2nd 2019 @ CNBC CMU

function plot_raw_raster(all_spike_train)


figure('Position', [100, 100, 1800, 800]);
for ii = 1:size(all_spike_train,1)
plot(find(all_spike_train(ii,:)), ii*ones(sum(all_spike_train(ii,:)),1),...
    's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'Markersize', 2); hold on
end
grid on; 







