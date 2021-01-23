function fitnessVal = evalFitness(x, bearing_diameter, ...
    sun_shaft_diameter, maxX, GR_front_target, GR_rear_target)

phi = fval(x, GR_front_target, GR_rear_target);
penalty = constraint(x, bearing_diameter, sun_shaft_diameter, maxX);

fitnessVal = phi + penalty;