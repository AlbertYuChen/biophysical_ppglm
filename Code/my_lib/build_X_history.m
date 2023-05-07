% Author: Yu Chen
% Date: Feb 23rd 2019 @ CNBC


function [Xsp, model] = build_X_history(all_spike_train, Trial_Index_train, model)
dt = model.dt;

% for ti = Trial_Index_train
%     spike_index_array = MC(neuron_index).spikeIndices == ti;
%     single_trail_spike_time = MC(neuron_index).spikeTimes(spike_index_array) * dt;
%     
%     HF_u = get_IMI_input_time(single_trail_spike_time, stim_time, 1);
%     zero_area = round(single_trail_spike_time/dt);
%     HF_u(1:zero_area) = NaN;      
%     HF_u = HF_u(valid_point_range);
%     HF_u(isnan(HF_u)) = [];
% end
% 
% HIST_bin_width = 0.001;  
% figure('Position', [300, 100, 600, 500]); 
% 
% histH = histogram(HF_u, 'Normalization', 'pdf', 'BinWidth', HIST_bin_width); 
% 
% ylabel('PDF', 'interpreter', 'latex')
% xlabel('time [sec]', 'interpreter', 'latex')
% 
% % ------- design H basis --------
% h_num_basis = 8;
% [BasisBegin, BasisEnd, NodeArry] = NodeArray_Distribution(histH, h_num_basis);
% B_s = 0;
% B_e = BasisEnd + 0.01;
% 
% spline_order = 4;
% knots_intense = [0 0 0 0 NodeArry B_e];
% % knots_intense = [0 0 0 0 0.0036 0.0076 0.0115 0.0154 0.0194 0.0241 0.0320 0.0890];
% % model.ht = (0:dt:0.089)';
% 
% model.ht = (B_s:dt:B_e)';
% 
% XBbas = [];
% for spline_index = 0 : numel(knots_intense) - spline_order - 1
%     [yy, x1] = bspline_basis(spline_index, spline_order, knots_intense, model.ht); 
%     XBbas = [XBbas yy];
% end
% 
% model.nh = size(XBbas, 2); 
% nh = model.nh;
% 
% model.hbas = XBbas;
% model.hbasExt = [XBbas; zeros(8000, model.nh)];
% model.knots_intense = knots_intense;

%-------------------------
nht = model.nht; % % ballah 2000, alon 250
model.ht = dt*(0:nht-1)';  % time relative to spike of stim filter taps
nh = model.nh;  % number of basis vectors for representing h

[~,model.hbas,~] = make_Pillow_Basis(nht, nh, 8);

%-------------------------
% figure('Position', [300, 300, 500, 400]);
% plot(model.ht, model.hbas)
% xlabel('time [sec]')
% title('Basis for post-spike history filter')

% ---- Create Spike-history Design Matrix ------------------------
Xsp = [];
for ti = Trial_Index_train
    Xsp_irt = hist_conv(all_spike_train(ti,:)', model.hbas);
    Xsp = [Xsp; Xsp_irt(model.valid_point_range, :)];
end


