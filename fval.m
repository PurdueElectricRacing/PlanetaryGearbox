
function val = fval(x, GR_front_target, GR_rear_target)
%% Gear Ratio Objective Function
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

%% Calculate gear ratios
front_gr = (front_planet_s1_teeth * front_ring_teeth) / ...
    (front_sun_teeth * front_planet_s2_teeth) + 1;

rear_gr = (rear_planet_s1_teeth * rear_ring_teeth) / ...
    (rear_sun_teeth * rear_planet_s2_teeth) + 1;


%% Calculate objective
val = (front_gr - GR_front_target) ^ 2;
val = val + (rear_gr - GR_rear_target) ^ 2;

