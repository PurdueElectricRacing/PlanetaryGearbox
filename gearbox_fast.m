function compoundgear
%% NOTE
%This function currently produces potential coumpound planetary gear
%setups for different diametral pitch values (defined below) and other
%input parameters. 
%
%% INPUTS
% input parameters
pd_stand = [16 18 20 24 28]; % diametral pitch
min_sun_diameter = 0.85;     % minimum pitch diameter of the sun gear
max_sun_diameter = 2;        % maximum pitch diameter of the sun gear
min_planet_diameter = 0.85;  % minimum pitch diameter of stage 2 planet gear
minOD = 2;                   % minimum outer diameter of gearbox, using pitch diameter of sun & stage 1 planet gear
maxOD = 4.5;                 % maximum outer diameter of gearbox, using pitch diameter of sun & stage 1 planet gear
planets = 3;                 % number of planet gears
gr = [6 6.5 7.5 8.5];        % range of gear ratios to be accepted for front and rear [front front rear rear]

%% CALCULATIONS
% generate initial matrices of all possible and somewhat realistic gears
% calculate max/min diameters for ring and planet gears
min_planet1_diameter = (minOD - max_sun_diameter) / 2;
max_planet1_diameter = (maxOD - min_sun_diameter) / 2;

min_planet2_diameter = min_planet1_diameter;
max_planet2_diameter = max_planet1_diameter;

min_ring_diameter = min_sun_diameter + (2 * min_planet2_diameter);
max_ring_diameter = maxOD;


% calculate approximate max/min number of teeth for each gear
min_sun_teeth = floor(pd_stand(1) * min_sun_diameter);
max_sun_teeth = ceil(pd_stand(end) * max_sun_diameter);

min_ring_teeth = floor(pd_stand(1) * min_ring_diameter);
max_ring_teeth = ceil(pd_stand(end) * max_ring_diameter);

min_planet1_teeth = floor(pd_stand(1) * min_planet1_diameter);
max_planet1_teeth = ceil(pd_stand(end) * max_planet1_diameter);

min_planet2_teeth = floor(pd_stand(1) * min_planet2_diameter);
max_planet2_teeth = ceil(pd_stand(end) * max_planet2_diameter);


% apply constraint that sun/ring gears teeth count must be a multiple of 3
while (min_sun_teeth / 3) ~= round(min_sun_teeth / 3)
   min_sun_teeth = min_sun_teeth - 1; 
end

while (max_sun_teeth / 3) ~= round(max_sun_teeth / 3)
   max_sun_teeth = max_sun_teeth + 1;
end

while (min_ring_teeth / 3) ~= round(min_ring_teeth / 3)
   min_ring_teeth = min_ring_teeth - 1; 
end

while (max_ring_teeth / 3) ~= round(max_ring_teeth / 3)
   max_ring_teeth = max_ring_teeth + 1; 
end


% calculate matrix of all possible gear diameters
sun_teeth = transpose([min_sun_teeth:3:max_sun_teeth]);
sun_diameters = sun_teeth ./ pd_stand;

planet1_teeth = transpose(min_planet1_teeth:max_planet1_teeth);
planet1_diameters = planet1_teeth ./ pd_stand;

ring_teeth = transpose([min_ring_teeth:3:max_ring_teeth]);
ring_diameters = ring_teeth ./ pd_stand;

planet2_teeth = transpose(min_planet2_teeth:max_planet2_teeth);
planet2_diameters = planet2_teeth ./ pd_stand;


% construct a matrix of all possible teeth configurations for the 4 gears
num_sun_gears = length(sun_teeth);
num_planet1_gears = length(planet1_teeth);
num_planet2_gears = length(planet2_teeth);
num_ring_gears = length(ring_teeth);

% matrix of valid setups
valid_setups = zeros(1,8);
num_pitches = length(pd_stand);
counter = 0;

% for loop iterating through every combination of gears. there are 3 checks
% to see if the setup if valid
for i = 1:num_pitches
    pd = pd_stand(i);
    for sun = 1:num_sun_gears
        for planet1 = 1:num_planet1_gears
            test_OD = sun_diameters(sun,i) + (2 * planet1_diameters(planet1,i));
            if (maxOD >= test_OD)
                for planet2 = 1:num_planet2_gears
                    for ring = 1:num_ring_gears
                        planet1_axis = sun_diameters(sun,i) + planet1_diameters(planet1,i);
                        planet2_axis = ring_diameters(ring,i) - planet2_diameters(planet2,i);
                        if planet1_axis == planet2_axis
                            gr_test = (planet1_teeth(planet1) / sun_teeth(sun) * ring_teeth(ring) / planet2_teeth(planet2)) + 1;
                            if (gr_test >= gr(1) && gr_test <= gr(2)) || (gr_test >= gr(3) && gr_test <= gr(4))
                                counter = counter + 1;
                                setup = [pd sun_teeth(sun) planet1_teeth(planet1) planet2_teeth(planet2) ring_teeth(ring) planet1_axis gr_test test_OD];
                                valid_setups = [valid_setups;setup];
                            end
                        end
                    end
                end
            end
        end
    end
end

% filter out planet 2 gear too small, plant 1 gear too big, sun gear too
% small
valid_setups = valid_setups(2:end,1:end);
i = 1;
length_1 = length(valid_setups);
counter = 1;

while counter <= length_1
    if valid_setups(i,3) >= (valid_setups(i,2) + valid_setups(i,3)) * sind(180/planets)
        valid_setups(i,:) = [];
    elseif (valid_setups(i,4) / valid_setups(i,1)) < min_planet_diameter
        valid_setups(i,:) = [];
    elseif (valid_setups(i,2) / valid_setups(i,1)) < min_sun_diameter
        valid_setups(i,:) = [];
    else
        i = i + 1;
    end
    counter = counter + 1;
end

% get rid of any path diameters that don't have a high and low gear ratio
elem = length(valid_setups);
valid_path_diameters = zeros(elem,1);
j = 1;

for i = 1:elem
   if valid_setups(i,7) >= gr(3)
       valid_path_diameters(j) = valid_setups(i,6);
       j = j + 1;
   end
end


i = 1;
while i <= length(valid_setups)
   if ismember(valid_setups(i,6), valid_path_diameters)
      i = i + 1; 
   else
       valid_setups(i,:) = [];
   end
end

xlswrite('gear_box.xlsx',valid_setups)


def = 0;
