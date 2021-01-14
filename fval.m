function val = fval(x, GR)
%% Planetary Gear Ratio Target Function
% Elliot Stockwell 1/12/21
%
% x : array of input variables (planet_d1, planet_d2, sun_d1, dp)
% GR : target gear ratio

planet_d1 = x(1);
planet_d2 = x(2);
sun_d1 = x(3);
ring_d2 = planet_d1 + planet_d2 + sun_d1;

val = (((planet_d1 * ring_d2 / sun_d1 / planet_d2) + 1) - GR)^2;