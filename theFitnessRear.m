function val = theFitnessRear(x, bearing_d, GR, planets_diameter_axis, ...
    sun_extrusion_d)
%% Planetary Gearbox GA Fitness Function
% Elliot Stockwell 1/12/21
%
% x : array of input variables (planet_d1, planet_d2, dp)
% bearing_d : inner diameter of the bearing
% GR : target gear ratio

y = zeros(4,1);
y(1) = x(1);
y(2) = x(2);
y(3) = planets_diameter_axis - y(1);
y(4) = x(3);

objective = fval(y, GR);
penalty = constraint(y, bearing_d, sun_extrusion_d);
val = objective + penalty;