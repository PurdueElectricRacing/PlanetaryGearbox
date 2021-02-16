function val = theFitnessFront(x, bearing_d, GR, sun_extrusion_d)
%% Planetary Gearbox GA Fitness Function
% Elliot Stockwell 1/12/21
%
% x : array of input variables (planet_d1, planet_d2, sun_d1)
% bearing_d : inner diameter of the bearing
% dp : diametric pitch
% GR : target gear ratio

objective = fval(x, GR);
penalty = constraint(x, bearing_d, sun_extrusion_d);
val = objective + penalty;