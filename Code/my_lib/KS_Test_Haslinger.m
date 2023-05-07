% Author: Yu Chen
% Date: May 19th 2018 @ CNBC

function Tau = KS_Test_Haslinger(Tau)

rng(623)
% rng('default')
zeta_list = [];

if ~isfield(Tau, 'zeta_list') && isfield(Tau, 'p_k_i')

for nn = 1:length(Tau.p_k_i)
    p_k_i_list = Tau.p_k_i{nn};
    q_k_i_list = -log( 1-p_k_i_list );
    zeta = sum( q_k_i_list(1:(end-1)) );
    r = rand(1);
    p_last = p_k_i_list(end);
    delta = -log( 1-r*p_last );
    zeta_list = [zeta_list, zeta+delta];
    
%     lambda_k_i = Tau.lambda_k_i{nn};
%     zeta = sum( lambda_k_i(1:(end-1)) );
%     r = rand(1);
%     p_last = 1-exp(-lambda_k_i(end));
%     delta = -log( 1 - r*p_last );
%     zeta_list = [zeta_list, zeta+delta];
end

Tau.zeta_list = zeta_list;
Tau.Uniform_list = 1- exp(-zeta_list);

end

zeta_list = Tau.zeta_list;

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
if ~isfield(Tau, 'U1')

Tau.nisi = Tau.nspk - 1;
endindex = cumsum(Tau.nisi);
startindex = endindex;
startindex = circshift(startindex, 1) + 1; 
startindex(1) = 1;

zeta_list1 = zeta_list;
zeta_list2 = zeta_list;
zeta_list1(endindex) = [];
zeta_list2(startindex) = [];

Ui1 = 1 - exp(-zeta_list1);
Ui2 = 1 - exp(-zeta_list2);

else
Ui1 = 1 - exp(-Tau.U1);
Ui2 = 1 - exp(-Tau.U2);
end

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












