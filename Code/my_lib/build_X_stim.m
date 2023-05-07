% Author: Yu Chen
% Date: Feb 23rd 2019 @ CNBC


function [Xstim, model] = build_X_stim(stim, Trial_Index_train, model)

dt = model.dt;

%-------------------------------
nkt = model.nkt; % ballah 80, alon 40
model.kt = dt*(0:nkt-1)';
nk = model.nk;

[~,model.kbas,~] = make_Pillow_Basis(nkt, nk, 10);

%-------------------------------
% Kbasis = zeros(10, 7);
% Kbasis(1,1) = 1;
% Kbasis(2,2) = 1;
% Kbasis(3,3) = 1;
% Kbasis(4,4) = 1;
% Kbasis(5:6,5) = 1;
% Kbasis(7:8,6) = 1;
% Kbasis(9:10,7) = 1;
% 
% nkt = 10;
% nk = 7;
% model.kt = dtStim*(0:nkt-1)';
% model.ktbas = Kbasis; 
%-------------------------------
% nkt = 80;
% spline_order = 4;
% % knots_intense = [0 0 0 0 0.001 0.002 0.003 0.005 0.007 0.01  0.015 0.03]; 
% knots_intense = [0 0 0   0:0.003:0.02   0.026:0.006:0.06   0.08 0.08];  % homo knots
% % knots_intense = [0 0 0 0 0.001 0.003 0.007 0.01  0.015 0.03]; 
% % knots_intense = [0 0 0 0 0.001 0.002 0.003 0.005 0.007 0.01  0.015 0.03 0.035 0.05]; 
% model.kt = (0:dt:nkt*dt)';
% 
% KBbas = [];
% for spline_index = 0 : numel(knots_intense) - spline_order - 1
%     [yy, x1] = bspline_basis(spline_index, spline_order, knots_intense, model.kt); 
%     KBbas = [KBbas yy];
% end
% model.ktbas = KBbas;
% nk = size(KBbas,2);

% -------------------------------
% figure('Position', [300, 300, 500, 400]);
% plot(model.kt, model.kbas)
% xlabel('time [sec]')
% title('Basis for stimulus filter')
% ---- Convolve stimulus with spatial and temporal bases -----
Xstim = [];
model.k_delay = 0;
for ti = Trial_Index_train
stim_cov_filterbasis = stim_conv(stim(:,ti), model.kbas, model.k_delay);
Xstim =  [Xstim; stim_cov_filterbasis(model.valid_point_range,:)];
end
























