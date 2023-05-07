% Author: Yu Chen
% Date: Feb 24th 2019 @ CNBC

function plot_models(model_list)

if length(model_list) == 1
model = model_list;

% ---------------------------- K -------------------------
% model.kw = prsML(1:nk);

figure('Position', [300, 300, 400, 250]);
plot(model.kt, model.k, 'r'); 
hold on;
xlabel('Time [sec]')
title('Stim filter');
% legend('K_0','K_1', 'location', 'northwest')
grid on
%------------------------- DC + KDC ----------------------

% y0 = stim_conv(stim, model.k, model.k_delay);
% 
% figure('Position', [300, 300, 400, 250]); 
% plot(stim, 'color', [.6 .6 .6]); hold on
% plot(y0, 'r')
% % load('D:\CRCNS_Project\Output\N1_Pillow_logit_1000-2000ms_22-Jun-2018')
% % y1 = stim_conv(stim, mod_set.k, mod_set.k_delay);
% % plot(y1, 'b')
% 
% xlim([600 900])
% xlabel('Time [ms]')
% % legend('Stim','Stim * K_0','Stim * K_1')
% title('Stim convolves with filter')
%---------------------------- H ---------------------------
% model.hw = prsML(nk+ndc+1:nk+nh+ndc);


figure('Position', [300, 300, 400, 250]);
plot(model.ht, model.h,'r'); hold on
% load('D:\CRCNS_Project\Output\N1_Pillow_logit_1000-2000ms_22-Jun-2018')
% plot(mod_set.ht, mod_set.h,'b');
title('Post-spike filter');
% legend('H_0','H_1', 'location', 'northwest')
ylim([-80 5])
xlim([0 model.ht(end)])
grid on

%---------------------- beta ------------------------------
figure
plot([model.kw; model.dc; model.hw], 'r.')

else % ++++++++++++++++++++++++++++++++++++++++++++

KW = [];
HW = [];
DC = [];

g = get_model_par(model_list, 'channel_scalar');
colorlist = jet( length(model_list) );
% ======================= K =======================

figure('Position', [300, 300, 500, 400]);
hold on;

for i = 1:length(model_list)

% x = gg0.ttk;
% y = gg0.k;
% dy = gg0.k_SE;
% fill([x;flipud(x)], [y-dy; flipud(y+dy)], [.9 .85 .85], 'linestyle','none');

plot(model_list{i}.kt, model_list{i}.k, 'color', colorlist(i,:), 'linewidth', 3);
KW = [KW, model_list{i}.kw];

end

% ylim([-2 6])

title('Stimulus filters');
ylabel('log firing rate')
xlabel('time [s]')
grid on
legend({num2str(g')})

% ======================= DC + KDC =======================

figure('Position', [300, 300, 500, 400]);
hold on;
for i = 1:length(model_list)

% x = gg0.iht;
% y = gg0.ih;
% dy = gg0.h_SE;
% fill([x;flipud(x)], [y-dy; flipud(y+dy)], [.9 .85 .85], 'linestyle','none');

plot(g(i), model_list{i}.dc, '+', 'color', colorlist(i,:), 'linewidth', 3);

DC = [DC, model_list{i}.dc];
end 

grid on
xlabel('channel conductance scalar')
% ylabel('beta value')
title('GLM parameter of baseline')
legend({num2str(g')})
% set(gca,'XTickLabel',{g})

% ======================= H =======================

figure('Position', [300, 300, 500, 400]);
hold on;

for i = 1:length(model_list)

% x = gg0.ttk;
% y = gg0.k;
% dy = gg0.k_SE;
% fill([x;flipud(x)], [y-dy; flipud(y+dy)], [.9 .85 .85], 'linestyle','none');

plot(model_list{i}.ht, model_list{i}.h, 'color', colorlist(i,:), 'linewidth', 3);

HW = [HW, model_list{i}.hw];

end



grid on
title('Post-spike filter');
xlabel('time [s]')
ylabel('log firing rate')
% ylim([-40 10])
xlim([0 model_list{i}.ht(end)])
% legend(legend_list, 'location', 'southeast')    
legend({num2str(g')})



%%
log_g = log(g);
x = g;

figure; hold on
plot(x, DC, '-o')
xlabel('channel conductance')
title('GLM parameter of baseline')

% ------------------- K -------------------
colorlist = jet( model_list{1}.nk );
figure; hold on
for i = 1:model_list{1}.nk
plot(x, KW(i,:), '-o', 'color', colorlist(i,:) )
end
grid on
xlabel('channel conductance')
title('GLM parameter of stimulus filter')


figure; hold on
for i = 1:model_list{1}.nk
plot(model_list{1}.kt, model_list{1}.kbas(:,i), 'color', colorlist(i,:) )
end
grid on
title('stimulus filter basis')
xlabel('time [s]')

% ------------------- H -------------------
colorlist = jet( model_list{1}.nh );
figure; hold on
for i = 1:model_list{1}.nh 
plot(x, HW(i,:), '-o', 'color', colorlist(i,:) )
end
grid on
xlabel('channel conductance')
title('GLM parameter of history filter')

figure; hold on
for i = 1:model_list{1}.nk
plot(model_list{1}.ht, model_list{1}.hbas(:,i), 'color', colorlist(i,:) )
end
grid on
title('history filter basis')
xlabel('time [s]')

end % end of the 












