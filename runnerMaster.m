%% Runs the runner function

%% Constants
pop_factor = 4;
num_front_designs = 5;
GR_front = 6.25;
GR_rear = 7.75;

%% Create Threadpool
p = gcp('nocreate');

if ~isempty(p)
    delete(p);
end

cpu_threads = max(1, str2double(getenv('NUMBER_OF_PROCESSORS')) - 2);
p = parpool(cpu_threads);

%% Run Front GA

num_runs = 100;
data_front = zeros(num_runs * num_front_designs, 11);

parfor i = 1:num_runs * num_front_designs
    [FOPT, planet_d1_act, planet_d2_act, sun_d1_act, ring_d2_act, ...
    dp, gr_act] = runnerFront(pop_factor, GR_front)
    data_front(i, :) = [FOPT, planet_d1_act, planet_d2_act, sun_d1_act, ...
                ring_d2_act, dp, gr_act, planet_d1_act * dp, ...
                planet_d2_act * dp, sun_d1_act * dp, ring_d2_act * dp];    
end

csvwrite('ga_results_front.csv', data_front);

[~, indexes] = sort(data_front(:, 1));
best_indexes = indexes(1:num_front_designs);
top_runs = data_front(best_indexes, :);

fprintf("\nFront GA Complete; Starting Rear GA\n")

%% Run Rear GA

data_rear = zeros(num_runs * num_front_designs, 11);
planets_diameter_axis_short = sum(top_runs(:, [2 4]), 2);
planets_diameter_axis_arr = zeros(num_runs * num_front_designs, 1);

for j = 1:num_runs * num_front_designs
    lookup_idx = floor((j - 1) / num_runs) + 1;
    planets_diameter_axis_arr(j) = planets_diameter_axis_short(lookup_idx);
end

parfor i = 1:num_runs * num_front_designs
    planets_diameter_axis = planets_diameter_axis_arr(i);
    
    [FOPT, planet_d1_act, planet_d2_act, sun_d1_act, ring_d2_act, ...
    dp, gr_act] = runnerRear(pop_factor, GR_rear, planets_diameter_axis)
    data_rear(i, :) = [FOPT, planet_d1_act, planet_d2_act, sun_d1_act, ...
                ring_d2_act, dp, gr_act, planet_d1_act * dp, ...
                planet_d2_act * dp, sun_d1_act * dp, ring_d2_act * dp];     
end

csvwrite('ga_results_rear.csv', data_rear);

%% Cleanup
delete(p)
