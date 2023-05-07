% Author: Yu Chen
% Date: Feb 23rd 2019 @ CNBC

function lamMax = get_lambdaMax(X_list, Y_list, beta, model_list)
B = length(model_list);
d = size(X_list{1}, 2);
n = size(X_list{1}, 1);
g_list = get_model_par(model_list, 'channel_scalar');
D = build_D(B, d, g_list);

for i = 1:B
mu = sigmoid(X_list{i} * beta(:,i));
d_beta(:,i) = X_list{i}' * (mu - Y_list{i}); 
end

x = (D*D')\D*d_beta(:);
lamMax = max(abs(x));

end


function D = build_D(B, d, g_list)
g = 1./diff(g_list);

D = sparse((B-1)*d, B*d);

for i = 1:(B-1)
    
    r = ((i-1)*d+1):i*d;
    c = r;
    matind = sub2ind([(B-1)*d, B*d], r, c);
    D(matind) = g(i);

    r = ((i-1)*d+1):i*d;
    c = r + d;
    matind = sub2ind([(B-1)*d, B*d], r, c);
    D(matind) = -g(i);

end

% figure
% imagesc(D)

end





























