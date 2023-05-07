% Author: Yu Chen
% Date: Aug 27th 2017 @ CNBC

function Tau = Tau_cross_validation(glmprs, stim, spktrMat, KS_valid_point_range, Itot_function)

N_trails = size(spktrMat, 1);

%% extract Tau, i.e. integral transformed intervals
Tau.Zscr = [];
Tau.Zi1 = [];
Tau.Zi2 = [];
Tau.wgt = [];

for ti = 1:N_trails
    spkInd = intersect(find(spktrMat(ti,:)), KS_valid_point_range);
    lambda = exp( Itot_function(glmprs, stim, spktrMat(ti,:)) );    
    Zscr_last = NaN;
    Tau.wgt = [Tau.wgt ones(1, spkInd(1))]; % add ones in the beginning of each trial. not the beginnig of the whole array.
    
    for ii=2:length(spkInd)
        Zscr_curr = sum(lambda((spkInd(ii-1)+1):spkInd(ii))); 
        Tau.Zscr = [Tau.Zscr Zscr_curr];
        Tau.wgt = [Tau.wgt  Zscr_curr*ones(1, spkInd(ii) - spkInd(ii-1))];
        
        %---- for iid test ----
        if ~isnan(Zscr_last)
        Tau.Zi1 = [Tau.Zi1 Zscr_last];
        Tau.Zi2 = [Tau.Zi2 Zscr_curr];
        end
        Zscr_last = Zscr_curr;
    end
    Tau.wgt = [Tau.wgt ones(1, length(KS_valid_point_range) - spkInd(end))];
end


















