function [Stress, Strength, FOS] = SpurGearBending(input1)
%% MOST CONSTANTS NOT DEFINED, STILL NEED TO SET UP FUNCTION INPUTS/OUTPUTS
% pitch circle
%
% diametral Pitch
% base circle
%% dP Equation %% 2 Independent Variables / 1 Independent Variables
% dP = pitch diameter pinion
% NP = Number of teeth on pinion
% Pd = Diametral pitch of pinion
dP = NP/Pd;
%% V Equation %% 2 Independent Variables / 1 Independent Variables
% V = 
% d = Operating pitch diameter of pinion
% n = Speed

V = pi*d*n/12;
%% Wt Equation %% 1 Independent Variables / 2 Dependent Variables
% Wt = Transmitted Load
% H = Power

Wt = 33000*H/V;
%% Kv Formula (Eq. 14-27, pg. 756) %% 2 Independent Variables / 3 Dependent Variables
% Kv = Dynamic Factor
% Qv = Transmission Accuracy-level number
% V = Magnitude of the pitch-line velocity?
B = 0.25*(12-Qv)^(2/3);
A = 50+56*(1-B);

Kv = ((A+sqrt(V))/A)^B; % for V in ft/min
% Kv = ((A+sqrt(200*V))/A)^B; % for V in m/s
%% Ks Formula (Eq.(a), section 14-30, pg. 759) % 1-3 Independent Variables / 1 Dependent Variable
% Ks = size factor (if Ks < 1 then use 1)
% P = Diametral Pitch
% F = Net Face width of narrowest member
Ks = 1; % or
% Ks = 1/kb or ?? what is kb?
% Ks = 1.192*(F*sqrt(Y)/P)^0.0535 ?? what is Y
%% Km Formula (Eq. 14-30, pg. 759) %% 7 Independent Variables / 4 Dependent Variables 
% F = Net Face Width of Narrowest Member
% Km = Load-Distribution Factor
if F <= 1 % inches
    Cpf = F/(10*d) - 0.025;
elseif F > 1 && F <= 17 
    Cpf = F/(10*d) - 0.0375 + 0.0125*F;
elseif F > 17 && F <= 40
    Cpf = F/(10*d) - 0.1109 + 0.0207*F - 0.000338*F^2;
end          
% Cmc = 1 (for uncrowned teeth) or  Cmc = 0.8 (for crowned teeth)
% Cpm = 1 (for straddle-mounted pinion w/ S1/S < 0.175) or = 1.1 (for straddle-mounted pinion w/ S1/S >= 0.175)
% A = ?, B = ?, C = ? see table 14-9 pg. 760
Cma = A + B*F + C*F^2;
% Ce = 0.8 (for gearing adjusted at assembly, compatibility improved by lapping, or both) or = 1 (for all other conditions)
Cmf = 1 + Cmc*(Cpf*Cpm+Cma*Ce); 
Km = Cmf;
%% KB Formula (Eq. 14-40, pg. 764) %% 3 Independent / 1 Dependent Variable
% KB = Rim-Thickness factor
% mB = backup ratio
% tR = rim thickness
% ht = tooth height

mB = tR/ht;

if mB < 1.2
    KB = 1.6*ln(2.242/mB);
else
    KB = 1; % mB >= 1.2
end
%% J Formula (Fig. 14-6, pg. 753) %% 4 Indepedent Variables / 1 Dependent
% J = Geometry factor for bending strength
% pN = normal base pitch
% Y = form factor
% Z = length of the line of action in the transverse plane
% Kf = fatigue stress-concentration factor
mN = pN/(0.95*Z); % load-sharing ratio
J = Y/(Kf*mN);
%% Stress Equation
% Stress = bending stress
% Ko = overload factor, see table on pg. 766
Stress = Wt*Ko*Kv*Ks*(Pd/F)*(Km*KB/J);
%% FOS %% ALL VARIABLES ARE INDEPENDENT
% St = AGMA bending strength (see table 14-3,14-4, pg. 748,749)
% YN = Stress cycle factor for bending strength (see Fig. 14-14, pg. 763)
% KT = Temperature Factor (1 if T < 250 degrees Fahrenheit
% KR = Reliability Factor (see Table 14-10, Eq. 14-38, pg. 763,764)
% FOS = Bending Factor of Safety

FOS = (St*YN/(KT*KR))/Stress;
%% Strength %% 4 Independent Variables / 1 Dependent
% St = AGMA bending strength (see table 14-3,14-4, pg. 748,749)
% SF = Safety Factor - bending
% YN = Stress cycle factor for bending strength (see Fig. 14-14, pg. 763)
% KT = Temperature Factor (1 if T < 250 degrees Fahrenheit)
% KR = Reliability Factor (see Table 14-10, Eq. 14-38, pg. 763,764)

Strength  = (St/FOS)*(YN/(KT*KR));
