% Author: Yu Chen
% Date: Feb 23rd 2019 @ CNBC

function par = get_model_par(model_cell_array, varName)

model_list = cell2mat(model_cell_array);

if isfield(model_list, varName)
    par = [model_list.(varName)];
else
   disp(['variable: ' varName ' does not exist']) 
end

