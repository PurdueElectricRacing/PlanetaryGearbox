function [FOPT, planet_d1_act, planet_d2_act, sun_d1_act, ring_d2_act, ...
    dp, gr_act] = runnerFront(pop_factor, GR)
%% Planetary Gearbox GA Runner Script
% Elliot Stockwell 1/12/21

%% Constants
adjustment_val = 0.2; %inches

planet_extrusion_d = 0.75 + adjustment_val; % inches
sun_extrusion_d = 0.71 + adjustment_val; % inches
bearing_d = 4.48; % inches

%% Executables
vlb = [planet_extrusion_d, planet_extrusion_d, sun_extrusion_d, 1];
vub = [0.45 * bearing_d, 0.45 * bearing_d, bearing_d, 4];
bits = [16, 16, 16, 2];
l = sum(bits);

parain = zeros(14, 1);
parain(1) = 0;
parain(11) = 4 * l * pop_factor;
parain(13) = (l + 1) / (2 * l * parain(11));
parain(14) = 200;

options = goptions(parain);

[X,FOPT,STATS,NFIT,FGEN,LGEN,LFIT] = GA550('theFitnessFront', ...
    [], options, vlb, vub, bits, bearing_d, GR, sun_extrusion_d);


%% Explain Results
planet_d1 = X(1);
planet_d2 = X(2);
sun_d1 = X(3);
dp = 14 + 2 * X(4);

%% Print Real Design
planet_d1_act = round(dp * planet_d1) / dp;
planet_d2_act = round(dp * planet_d2) / dp;
sun_d1_act = round(dp * sun_d1) / dp;
ring_d2_act = planet_d1_act + planet_d2_act + sun_d1_act;
gr_act = (planet_d1_act * ring_d2_act / sun_d1_act / planet_d2_act) + 1;

end