% Author: Yu Chen
% Date: Oct-06-2018 @ MI

function hist_Z_nind(Tau)

interval_Trial_Start_End = cell2mat(Tau.interval_Trial_Start_End');
normalized_index_list = interval_Trial_Start_End(:,4) ./ interval_Trial_Start_End(:,5); 
Uniform_list = Tau.Uniform_list;


% dim 1 (Y axis), dim 2 (X axis)
X = [Uniform_list', normalized_index_list];
% [N,c] = hist3(X, 'nbins', [10, 10]);
[N,c] = hist3(X, 'ctrs', {0.05:0.1:0.95, 0.05:0.1:0.95});

N = N / length(X) * 100;

figure
imagesc(N, [0, 3]);
set(gca, 'XTick', 1:10, 'XTickLabel', c{2})
set(gca, 'YTick', 1:10, 'YTickLabel', c{1})
xlabel('Normalized index')
ylabel('Z')
colorbar;


figure
imagesc(N);
set(gca, 'XTick', 1:10, 'XTickLabel', c{2})
set(gca, 'YTick', 1:10, 'YTickLabel', c{1})
xlabel('Normalized index')
ylabel('Z')
colorbar;



R = corrcoef(Uniform_list', normalized_index_list)


% figure
% hist3(X, 'CdataMode','auto')
% view(2)
% 
% figure 
% hist3(X, 'ctrs', {0.05:0.1:0.95, 0.5:1:9.5}, 'CdataMode','auto') 
% view(2)








