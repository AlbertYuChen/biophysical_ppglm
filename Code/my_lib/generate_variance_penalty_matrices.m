% Author: Yu Chen
% Date: Feb 26th 2017 @ CNBC
% The Omiga is used for smoothing spline regression

function Omiga = generate_variance_penalty_matrices(varargin)

% default values
normal = 1;
baseline = 0;

knot_blk_array = {};
ind = 1;

while ind <= nargin
    if ischar(varargin{ind})
        switch lower(varargin{ind})
        case 'normal'
            normal = varargin{ind+1};
            ind = ind + 2;
            continue;
        case 'baseline'
            baseline = varargin{ind+1};
            ind = ind + 2;
            continue;
        end      
    end
    
    if isvector(varargin{ind})
        knot_blk_array = {knot_blk_array varargin{ind}};
    end
    ind = ind + 1;
end

spline_order = 4;  % 4 is cubic, 3 is quadratic

if baseline > 0
    Omiga = baseline;
else
    Omiga = [];
end
for blk_ind = 2:length(knot_blk_array)
    
    knot_list = knot_blk_array{blk_ind};
    basis_time = knot_list(1):0.001:knot_list(end);    
%     figure;
%     hold on;   
    spline_basis = [];
    for spline_index = 0 : numel(knot_blk_array{blk_ind}) - spline_order - 1
        [yy, x1] = bspline_basis(spline_index, spline_order, knot_blk_array{blk_ind}, basis_time); 
        spline_basis = [spline_basis yy'];
%         plot(x1, yy, '.');
    end
    
    spline_basis = diff(diff(spline_basis));
    if normal
        spline_basis = normc(spline_basis);
    end
    Omiga = blkdiag(Omiga, spline_basis'*spline_basis);
end

% figure
% colormap('hot'); 
% imagesc(Omiga);
