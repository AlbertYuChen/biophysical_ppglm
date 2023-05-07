% Author: Yu Chen
% Date: Jul 25th 2018 @ Beijing


function out = sigmoid(x, clip)

if nargin<2
    clip = true;
end

if ~clip
    out = 1 ./ (1 + exp(-x));
else
    dataClass = class(x);
    lowerBnd = log(eps(dataClass)); 
    upperBnd = -lowerBnd;
    out = 1 ./ (1 + exp(-constrain(x,lowerBnd,upperBnd)));
end

end

function x = constrain(x,lower,upper)
% Constrain between upper and lower limits, and do not ignore NaN
x(x<lower) = lower;
x(x>upper) = upper;
end