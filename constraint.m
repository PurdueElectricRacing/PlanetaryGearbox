function val = constraint(x, bearing_d, sun_extrusion_d)
%% Planetary Gearbox GA Constraint Function
% Elliot Stockwell 1/12/21
%
% x : array of input variables (planet_d1, planet_d2, sun_d1, dp)
% bearing_d : inner diameter of the bearing

pm = 1E2; % penalty multiplier

planet_d1 = x(1);
planet_d2 = x(2);
sun_d1 = x(3);
dp = 12 + 4 * x(4);
ring_d2 = planet_d1 + planet_d2 + sun_d1;

g = zeros(1, 8);

mod_tolerance = .01;

%% Make sure teeth are integers
planet_1_mod = mod(dp * planet_d1, 1);
planet_2_mod = mod(dp * planet_d2, 1);
sun_mod = mod(dp * sun_d1, 3);
ring_mod = mod(dp * ring_d2, 3);

if sun_mod < 3 - mod_tolerance && sun_mod > mod_tolerance
    g(1) = pm;
end

if ring_mod < 3 - mod_tolerance && ring_mod > mod_tolerance
    g(2) = pm;
end

if planet_1_mod < 1 - mod_tolerance && planet_1_mod > mod_tolerance
    g(3) = pm;
end

if planet_2_mod < 1 - mod_tolerance && planet_2_mod > mod_tolerance
    g(4) = pm;
end

%% Side constraints
if dp > 24
    g(5) = pm;
end

g(6) = max(0, (sun_d1 + 2 * planet_d1) / bearing_d - 1) * pm;

g(7) = max(0, (planet_d2 / (0.45 * ring_d2)) - 1) * pm;

if sun_d1 < sun_extrusion_d
    g(8) = pm;
end
%% Sum constriants
val = sum(g);


