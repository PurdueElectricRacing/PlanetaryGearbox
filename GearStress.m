function GearStress(pd, ring_n, sun_n, planet_n, J_s, J_p)
%% NOTE
%INPUTS:
%pd = diametral pitch (don't confuse this as pitch diameter. Following
%notation from Shigleys)
%ring_n = number of teeth on ring gear
%sun_n = number of teeth on sun gear
%planet_n = number of teeth on planet gear
%
%J_s and J_p:
%Due to complexity of curve, the geometry factor (J) is determined outside
%of the code using Figure 14-6 from Shigleys 10th ed
%
%This code does not check for valid inputs of run, sun, and planet teeth.
%These must be figured out using the GearCalculator.m script

%% PREDEFINED CONSTANT
%Track and Vehicle Variables
end_time = 25; %length of endurance in minutes
Ne = 160; %Number of endurance events ran
motor_rpm = 4000; %[rpm]
torque = 15; %[Nm]
hp = torque*motor_rpm/9.5488/1000*1.34; %[hp]

%Gear Variables
dp_s = sun_n/pd; %pitch diameter sun [in]
dp_p = planet_n/pd; %pitch diameter planet [in]
Pa = 20; %Pressure angle [deg]
F_s = 1; %face width of sun [in]
F_p = 1; %face width of planet [in]

Y_tab = xlsread('Lewis Factor'); %Lewis Form Factor Table
N_s = motor_rpm * end_time * Ne * 3; %load cycle for sun gear
N_p = motor_rpm / (planet_n/sun_n) * end_time * Ne; %load cycle for planet gear

%Material/Manufacturing Properties
qV = 12; %quality 
% v_s %poisson's ratio for sun
% v_p %poisson's ratio for planet
% E_s %Young's modulus of sun
% E_p %Young's modulus of planet
Hb = 400; %Brinell Hardness guess

%% STRESS COEFFICIENT CALCULATIONS
V = (pi*dp_s*motor_rpm)/12; %Pitch Line Velocity
Wt = 33000*hp/V; %Transmitted Load

%Kv - Dynamic Factor
B = 0.25*(12-qV)^(2/3);
A = 50 + 56*(1-B);
Kv = ((A + sqrt(V))/A)^B;

%Ks - Size Factor
teeth = [sun_n planet_n];
Y = [0 0];
for j=1:2
    for i=1:26
        if Y_tab(i,1) == teeth(j)
            Y(j) = Y_tab(i,2);
            break
        elseif teeth(j) > Y_tab(i,1) && teeth(j) < Y_tab(i+1,1)
            Y(j) = Y_tab(i,2);
            break
        end
    end
end
Ks_s = 1.192*(F_s*sqrt(Y(1))/pd)^0.0535;
Ks_p = 1.192*(F_p*sqrt(Y(2))/pd)^0.0535;

%Load Distribution Factor
Cmc = 1; %assume uncrowned (.8 if crowned)
if F_s <= 1
    Cpf = F_s/10/dp_s-.025;
else
    Cpf = F_s/10/dp_s - 0.0375 + 0.0125*F_s;
end
Cpm = 1; %assume S1/s<.175
Cma= .00360 + .0102*F_s + (-.822*10^-4)*F_s^2; %assume precision
Ce = 0.9; %other conditions
Km = 1 + Cmc * (Cpf*Cpm + Cma*Ce);

%Rim Thickness Factor
Kb=1; %assume thick enough rim

%Overload Factor
Ko = 1.5;

%Bending Stress Cycle Factor
if N_s < 10^7
    Yn_s = 3.517*N_s^(-0.0617);
else
    Yn_s = 1.3558*N_s^(-0.0178);
end

if N_p < 10^7
    Yn_p = 3.517*N_p^(-0.0617);
else
    Yn_p = 1.3558*N_p^(-0.0178);
end

%Stress-cycle Factor Equations
if N_s < 10^7
    Zn_s = 1.249*N_s^(-0.0138);
else
    Zn_s = 1.4488*N_s^(-0.023);
end

if N_p < 10^7
    Zn_p = 1.249*N_p^(-0.0138);
else
    Zn_p = 1.4488*N_p^(-0.023);
end

%Reliability Factor
Kr = 1; %.99 Reliability

%Temperature Factor 
Kt = 1;

%Surface Condition Factor
Cf = 1; %Subject to change from manufacturer

%Load Sharing Ratio
mN = 1;

%Pitting Resistance Geometry Factor
mG = planet_n/sun_n;
I = (cosd(Pa)*sind(Pa))/(2*mN)*mG/(mG+1);

%Elastic coefficient
%this value can be either found using Eq 14-13 or Table 14-8
Cp = 2300; %Steel to steel [psi^(1/2)]
%Cp = (pi*((1-v_s^2)/E_s + (1-v_p^2)/E_p))^(-0.5);

%AGMA Strength Equations (Assume material same for sun and planet rn)
%Uncomment either grade 1 or grade 2 equations
%Grade 1
%St = 77.3*Hb + 12800; %Allowable bending stress number [psi]
%Sc = 322*Hb + 29000; %Allowable contact stress number [psi]
%Grade 2
St = 102*Hb + 16400; %Allowable bending stress number [psi]
Sc = 349*Hb + 34300; %Allowable contact stress number [psi]

%Hardness Ratio Factor
Ch = 1;

%% STRESS AND ENDURANCE CALCULATIONS
%%Only sun is being analzyed right now as sun is to fail before planet

%%Gear Wear
sigmaC_s = Cp*sqrt(Wt*Ko*Kv*Ks_s*Km/(dp_s*F_s)*Cf/I);
sH_s = Sc*Zn_s/(Kt*Kr)/sigmaC_s;

%Bending
sigma = Wt*Ko*Kv*Ks_s*pd/F_s*(Km*Kb)/J_s;
sF_s = St*Yn_s/(Kt*Kr)/sigma;

fprintf('Sun ring has max gear contact stress = %f psi and a FOS = %f\n',sigmaC_s, sH_s)
fprintf('Sun ring has max gear bending stress = %f psi and a FOS = %f\n',sigma, sF_s)

end