%% Master running file for GA methods

%% Define constants
bearing_diameter = 4.50;        % in
sun_shaft_diameter = 1.23;      % in
planet_shaft_diameter = 0.85;   % in

GR_front_target = 6.5;      % dimensionless
GR_rear_target = 8.00;      % dimensionless

numGenerations = 200;
num_runs = 500;

%% Define bounds
%
% X is formatted
% [ 
%   front_ring_teeth / 3
%   front_planet_s2_teeth
%   front_sun_teeth / 3
%   rear_planet_s2_teeth
%   rear_sun_teeth / 3
%   diametral_pitch
% ]

% Remember some of these are / 3
min_x = [...
    14, ...
    14, ...
    5, ...
    14, ...
    5, ...
    1];

max_x = [...
    35, ...
    48, ...
    35, ...
    48, ...
    35, ...
    4];


%% Perform last calculations
vlb = min_x;

vdiff = max_x - min_x;
bits = ceil(log(vdiff) / log(2));

vub = vlb + 2 .^ bits - 1;
maxX = max_x;


%% Run GA
data = zeros(num_runs, 20);

fprintf('\nBeginning GA\n');

parfor i=1:num_runs
    [FOPT, X] = runOneOptimization(vlb, vub, bits, bearing_diameter, ...
        planet_shaft_diameter, sun_shaft_diameter, maxX, ...
        GR_front_target, GR_rear_target, numGenerations)

    [front_diameters, front_teeth, rear_diameters, rear_teeth, ...
        dp, frontGR, rearGR] = explainX(X);
    
    data(i, :) = [FOPT, front_diameters, front_teeth, ...
        rear_diameters, rear_teeth, dp, frontGR, rearGR];
end


%% Save Info
filename = 'ga_results.csv';

header = {'Badness', ...
    'Front Planet S1 Diameter (in)', 'Front Planet S2 Diameter (in)', ...
    'Front Sun Diameter (in)', 'Front Ring Diameter (in)', ...
    'Front Planet S1 Teeth', 'Front Planet S2 Teeth', ...
    'Front Sun Teeth', 'Front Ring Teeth', ...
    'Rear Planet S1 Diameter (in)', 'Rear Planet S2 Diameter (in)', ...
    'Rear Sun Diameter (in)', 'Rear Ring Diameter (in)', ...
    'Rear Planet S1 Teeth', 'Rear Planet S2 Teeth', ...
    'Rear Sun Teeth', 'Rear Ring Teeth', ...
    'Diametral Pitch (1/in)', 'Front GR', 'Rear GR'};
    
commaHeader = [header;repmat({','},1,numel(header))];
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader);

fid = fopen(filename, 'w');
fprintf(fid, '%s\n', textHeader(1:end-1));
fclose(fid);

data_sorted = sortrows(data);
dlmwrite(filename, data_sorted, '-append');


%% Clean up
fprintf('\nDone.\n');

%% Create histogram

feasible_region_front = [6.25 6.25 6.75 6.75 6.25];
feasible_region_rear = [7.875 8.5 8.5 7.875 7.875];

figure(1);
hold on;grid on;grid minor;
colormap jet;
set(gcf, 'Position', [100 100 900 700]);
sct = scatter(data(:, end-1), data(:, end), 25, data(:, 1));
trg = plot(GR_front_target, GR_rear_target);
region = plot(feasible_region_front, feasible_region_rear);
set(trg, 'Marker', '.', 'MarkerSize', 25, 'Color', 'r', ...
    'LineStyle', 'None');
set(region, 'Color', 'k', 'LineStyle', '-.', 'LineWidth', 1);
title('Distribution of Results');
xlabel('Front Gear Ratio');ylabel('Rear Gear Ratio');
legend([sct, trg], 'GA Results', 'Target');


