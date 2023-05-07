% Author: Yu Chen
% Date: April 22nd 2017 @ CNBC
% Non-linear convolution using Volterra formular
% all input has to be column vectors

% g = [2 5 6]';
% out = x(t-1)*x(t-4)*x(t-5)


function out = convVolterra(x, g)
g = g(g~=0);
x = reshape(x,[],1);
g = reshape(g,[],1);

winlen = max(g);
xlen = length(x);
x = [zeros(winlen - 1, 1); x];
out = zeros(xlen, 1);
winsld = 1:winlen;

for ii = 1:xlen
    winval = x(winsld + ii - 1);
    winval = flipud(winval);
    out(ii) = prod(winval(g));
end
