function [FOPT, X] = runOneOptimization(vlb, vub, bits, bearing_diameter, ...
    planet_shaft_diameter, sun_shaft_diameter, maxX, GR_front_target, ...
    GR_rear_target, numGenerations)

% Chromosome length
l = sum(bits);

%% Define GA options
parain = zeros(14, 1);
parain(1) = 0;
parain(11) = 4 * l;
parain(13) = (l + 1) / (2 * l * parain(11));
parain(14) = numGenerations;

options = goptions(parain);

[X, FOPT, ~, ~, ~, ~, ~] = GA550('evalFitness', ...
    [], options, vlb, vub, bits, bearing_diameter, ...
    planet_shaft_diameter, sun_shaft_diameter, maxX, GR_front_target, ...
    GR_rear_target);