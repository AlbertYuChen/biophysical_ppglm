% Author: Yu Chen
% Date: May 15th 2018 @ CNBC

function score = PSTH_goodness_of_fit(PSTH0, PSTH_in, type)

if nargin < 3
    type = 'MSE';
end

switch type
    case 'MSE'
        score = mean((PSTH_in - PSTH0).^2);
    case 'KL'
        disp('KL')
        score = KLDiv(PSTH0, PSTH_in);
    otherwise
        disp('Wrong type')
end





