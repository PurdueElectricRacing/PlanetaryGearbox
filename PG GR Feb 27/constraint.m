function val = constraint(x, bearing_diameter, planet_shaft_diameter, sun_shaft_diameter, maxX)
%% Gear Ratio Constraint Function
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

%% GA Constants
PENALTY = 1e2;
val = 0;

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


%% Calculate side constraints
% Variables
comp_planet_diameter = (front_sun_diameter + front_planet_s1_diameter) / (rear_sun_diameter + rear_planet_s1_diameter) - 1;

% Don't exceed bearing bounds
if 2 * front_planet_s1_diameter + front_sun_diameter > bearing_diameter
    val = val + PENALTY;
end

if 2 * rear_planet_s1_diameter + rear_sun_diameter > bearing_diameter
    val = val + PENALTY;
end

if comp_planet_diameter < -0.01 || comp_planet_diameter > 0.01
    val = val + PENALTY;
end

% Make sure every gear is larger than the shaft that runs through it
check_planet_arr = [...
    front_planet_s1_diameter, ...
    front_planet_s2_diameter, ...
    rear_planet_s1_diameter, ...
    rear_planet_s2_diameter];

check_sun_arr = [...
    front_sun_diameter, ...
    rear_sun_diameter];

% Count occurances where the constriant is violated
val = val + PENALTY * sum(check_planet_arr < planet_shaft_diameter);

val = val + PENALTY * sum(check_sun_arr < sun_shaft_diameter);

% Enforce upper bounds on x
val = val + PENALTY * sum(x > maxX);
