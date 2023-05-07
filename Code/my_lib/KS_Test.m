function KS_Test(Tau)

[eCDF, zvals] = ecdf(Tau.Zscr); 
zvals(1)=0;
mCDF = 1-exp(-zvals); %Model CDF at z values.

%------ KS -------
figure('Position', [200, 200, 1400, 300])
subplot(141)
plot(mCDF, eCDF, 'Linewidth', 3)
hold on								%Freeze graphics.

% 0.95% CI
plot([0 1], [0 1]+1.36/sqrt(length(Tau.Zscr)), '--', 'Color', [.7 .7 .7])		%Upper confidence bound.
plot([0 1], [0 1]-1.36/sqrt(length(Tau.Zscr)), '--', 'Color', [.7 .7 .7])		%Lower confidence bound.

% 0.999% CI
% plot([0 1], [0 1]+1.95/sqrt(length(Tau.Zscr)), '--', 'Color', [.7 .7 .7])		%Upper confidence bound.
% plot([0 1], [0 1]-1.95/sqrt(length(Tau.Zscr)), '--', 'Color', [.7 .7 .7])		%Lower confidence bound.

hold off							%Release graphics window.
xlabel('Model CDF')					%Label the axes.
ylabel('Empirical CDF')
axis([0 1 0 1]);
title('KS test')
grid on; set(gca,'FontSize', 12)

% KS test p-value
epsilon95 = 1.36/sqrt(length(Tau.Zscr));
KS_p_value95 = 2*exp(-2*length(Tau.Zscr)*epsilon95^2);
epsilon = max(mCDF - eCDF);
KS_p_value = 2*exp(-2*length(Tau.Zscr)*epsilon^2);

disp(['epsilon  95%: ' num2str(epsilon95) '  p-value: ' num2str(KS_p_value95)])
disp(['epsilon data: ' num2str(epsilon) '  p-value: ' num2str(KS_p_value)])

%--- iid test ---
Ui1 = 1 - exp(-Tau.Zi1);
Ui2 = 1 - exp(-Tau.Zi2);
subplot(142)
plot(Ui1, Ui2, '.', 'Markersize', 4)
xlabel('U_i')
ylabel('U_{i+1}')
title('IID test')
axis([0 1 0 1])
set(gca,'FontSize', 12)

%------ histogram -------
subplot(143)
t = 0:0.01:10;
exppdf = exp(-t);
plot(t, exppdf, 'Color', [0 0.7 0], 'Linewidth', 2)
hold on
histogram(Tau.Zscr, 'Normalization', 'pdf', 'FaceColor', [0 0.4470 0.7410], 'BinWidth', 0.1)
xlim([0 5])
legend('Ideal curve')
title('PDF of \tau')
grid on; set(gca, 'FontSize',12)

subplot(144)
plot([0 1],[1 1], 'Color', [0 0.7 0], 'Linewidth', 2)
hold on
histogram(1- exp(-Tau.Zscr), 'Normalization', 'pdf', 'FaceColor', [0 0.4470 0.7410], 'BinWidth', 0.02)
xlim([0 1])
legend('Ideal curve', 'location', 'southeast')
title('PDF of 1 - exp\{-\tau\}')
grid on; set(gca, 'FontSize',12)

% figure 
% t = 0:0.01:10;
% exppdf = exp(-t);
% plot(t, exppdf, 'Color', [0 0.7 0], 'Linewidth', 2)
% hold on
% histogram(Tau.Zscr, 'Normalization', 'pdf', 'BinWidth', 0.05)
% xlim([0 4])
% % ylim([0 1])
% legend('Ideal curve')
% title('PDF of \tau')
% grid on; set(gca, 'FontSize',12)

%% for JNC paper usage
% figure('Position', [300, 100, 600, 500]);
% plot(mCDF, eCDF, 'Linewidth', 3)
% hold on								%Freeze graphics.
% 
% % 0.95% CI
% plot([0 1], [0 1]+1.36/sqrt(length(Tau.Zscr)), '--', 'Color', [.7 .7 .7]) 
% h_l = plot([0 1], [0 1]-1.36/sqrt(length(Tau.Zscr)), '--', 'Color', [.7 .7 .7]); 
% 
% xlabel('Model CDF', 'interpreter', 'latex') 
% ylabel('Empirical CDF', 'interpreter', 'latex')
% legend(h_l, {'95\% CI'}, 'Interpreter','latex', 'location', 'southeast');
% 
% axis([0 1 0 1]);
% % title('KS test', 'interpreter', 'latex')
% % title('\textit{Monkey-PMv} KS test', 'interpreter', 'latex')
% % title('\textit{Human-Cortex} homogeneous IMI KS test', 'interpreter', 'latex')
% set(gca, 'FontSize',18, 'TickLabelInterpreter', 'latex');

%% group activities
% figure('Position', [300, 100, 600, 500]); 
% load(['D:\Divergent_Spiketrain\Output\Monkey_LimitSpk_out.mat'])
% Tau = mod_set.Tau;
% [eCDF, zvals] = ecdf(Tau.Zscr);
% mCDF = 1-exp(-zvals);				%Model CDF at z values.
% h1 = plot(mCDF, eCDF, 'Linewidth', 3);
% hold on 
% 
% load(['D:\Divergent_Spiketrain\Output\Monkey_SpkRect_out.mat'])
% Tau = mod_set.Tau;
% [eCDF, zvals] = ecdf(Tau.Zscr);
% mCDF = 1-exp(-zvals);				%Model CDF at z values.
% h2 = plot(mCDF, eCDF, 'Linewidth', 3);
% hold on 

% 
% load(['D:\Divergent_Spiketrain\Output\Monkey_IMI_3h_out.mat'])
% Tau = mod_set.Tau;
% [eCDF, zvals] = ecdf(Tau.Zscr);
% mCDF = 1-exp(-zvals);				%Model CDF at z values.
% h3 = plot(mCDF, eCDF, 'Linewidth', 3, 'color', [0.4660    0.6740    0.1880]);
% hold on 
% 

% % 0.95% CI
% plot([0 1], [0 1]+1.36/sqrt(length(Tau.Zscr)), '--', 'Color', [.7 .7 .7]) 
% h_ci = plot([0 1], [0 1]-1.36/sqrt(length(Tau.Zscr)), '--', 'Color', [.7 .7 .7]); 
% 
% xlabel('Model CDF', 'interpreter', 'latex') 
% ylabel('Empirical CDF', 'interpreter', 'latex')
% legend([h1 h2 h3 h_ci], {'IMI-1','IMI-2','IMI-3','95\% CI'}, 'Interpreter','latex', 'location', 'southeast');
% 
% axis([0 1 0 1]);
% title('\textit{Monkey-PMv} IMI-1,2,3 models KS test', 'interpreter', 'latex')
% set(gca, 'FontSize',18, 'TickLabelInterpreter', 'latex');






















