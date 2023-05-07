
function u = get_time_stairs(cat_trail_spike_time_stairs, order)

N_trails = length(cat_trail_spike_time_stairs) / 2000;
u = zeros(1, length(cat_trail_spike_time_stairs));

for stim_section = 0:(N_trails - 1)
    max_ary = zeros(order, 1); 
    u(stim_section*2000 + 1) = 0.001;
    for bin_ind = 2:2000    
        if cat_trail_spike_time_stairs(stim_section*2000 + bin_ind - 1) > max_ary(1) + 0.0000001   
            max_ary = circshift(max_ary, 1);
            max_ary(1) = cat_trail_spike_time_stairs(stim_section*2000 + bin_ind - 1);   
            u(stim_section*2000 + bin_ind) = bin_ind / 1000  - max_ary(order);            
        else   
            u(stim_section*2000 + bin_ind) = bin_ind / 1000  - max_ary(order);           
        end
    end
end
