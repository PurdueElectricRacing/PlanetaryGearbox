function [front_diameters, front_teeth, ...
    rear_diameters, rear_teeth, dp, frontGR, rearGR] = explainX(x)

%% This function explains the design represented by x

% X is formatted
% [ 
%   front_ring_teeth / 3
%   front_planet_s2_teeth
%   front_sun_teeth / 3
%   rear_planet_s2_teeth
%   rear_sun_teeth / 3
%   diametral_pitch
% ]

%% Interpret x
front_ring_teeth = 3 * x(1);
front_planet_s2_teeth = x(2);
front_sun_teeth = 3 * x(3);
rear_planet_s2_teeth = x(4);
rear_sun_teeth = 3 * x(5);
dp = 12 + 4 * x(6);


%% Calculate derived variables
rear_ring_teeth = front_ring_teeth;

front_planet_s1_teeth = front_ring_teeth - ...
    front_planet_s2_teeth - front_sun_teeth;

rear_planet_s1_teeth = rear_ring_teeth - ...
    rear_planet_s2_teeth - rear_sun_teeth;


%% Calculate diameters
front_ring_diameter = front_ring_teeth / dp;
front_planet_s1_diameter = front_planet_s1_teeth / dp;
front_planet_s2_diameter = front_planet_s2_teeth / dp;
front_sun_diameter = front_sun_teeth / dp;

rear_ring_diameter = rear_ring_teeth / dp;
rear_planet_s1_diameter = rear_planet_s1_teeth / dp;
rear_planet_s2_diameter = rear_planet_s2_teeth / dp;
rear_sun_diameter = rear_sun_teeth / dp;


%% Prepare outputs
front_diameters = [front_planet_s1_diameter, front_planet_s2_diameter, ...
    front_sun_diameter, front_ring_diameter];

rear_diameters = [rear_planet_s1_diameter, rear_planet_s2_diameter, ...
    rear_sun_diameter, rear_ring_diameter];

front_teeth = [front_planet_s1_teeth, front_planet_s2_teeth, ...
    front_sun_teeth, front_ring_teeth];

rear_teeth = [rear_planet_s1_teeth, rear_planet_s2_teeth, ...
    rear_sun_teeth, rear_ring_teeth];

frontGR = (front_planet_s1_teeth * front_ring_teeth) / ...
    (front_sun_teeth * front_planet_s2_teeth) + 1;

rearGR = (rear_planet_s1_teeth * rear_ring_teeth) / ...
    (rear_sun_teeth * rear_planet_s2_teeth) + 1;


