function [Itot_NLN_Fold, SegInd, HeadsInd, NadirsInd] = findLocalNadirs(Itot_NLN,Itot,  N_trails, cutThreshold)

SegInd = {};
HeadsInd = {};
NadirsInd = {};

Itot_NLN_Fold = reshape(Itot_NLN, [], N_trails);
Itot_Fold = reshape(Itot, [], N_trails);

for ii = 1:N_trails
% for ii = 1
SegInd{ii} = find(Itot_NLN_Fold(:,ii) < cutThreshold);
IndTmp = [SegInd{ii}; 2002];
ContinuSeg = [0; diff(IndTmp)];
ContinuSegHead = find(ContinuSeg ~= 1);
NadirsInd{ii} = [];

    for jj = 1:(length(ContinuSegHead)-1)
        segRang = IndTmp( ContinuSegHead(jj):(ContinuSegHead(jj+1)-1) );
        [~, minInd] = min(Itot_NLN_Fold(segRang,ii));
        NadirsInd{ii} = [NadirsInd{ii}; segRang(minInd)];
    end
HeadsInd{ii} = IndTmp( ContinuSegHead(1:(end-1)) );
end


TrialInd = 1;

figure('Position', [300, 300, 1500, 350]);
plot(Itot_Fold(:, TrialInd), 'b');
hold on
plot(Itot_NLN_Fold(:, TrialInd), 'r');
% plot(SegInd{TrialInd}, Itot_NLN_Fold(SegInd{TrialInd}, TrialInd), 'co')
% plot(HeadsInd{TrialInd}, Itot_NLN_Fold(HeadsInd{TrialInd}, TrialInd), 'c+')
plot(NadirsInd{TrialInd}, Itot_NLN_Fold(NadirsInd{TrialInd}, TrialInd), 'gx')
% title('f_{Nonlinear}(X \ast b)'); 
% xlim([1 600])
grid on
ylim([-25 2])
legend('X \ast b', 'f_{Nonlinear}(X \ast b)')
set(gca, 'FontSize', 12)

HIST_bin_width = 0.5;
figure('Position', [200, 200, 1200, 300]);
histH_NLN = histogram(Itot_NLN, 100, 'Normalization', 'pdf');

figure('Position', [300, 300, 1500, 350]);
plot(Itot_Fold(:, TrialInd), 'b');
figure('Position', [300, 300, 1500, 350]);
plot(Itot_NLN_Fold(:, TrialInd), 'r');

