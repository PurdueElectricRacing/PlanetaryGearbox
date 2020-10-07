function [output1, output2] = SpurGearBending(input1)

%%
% dP = pitch diameter pinion
% NP = Number of teeth on pinion
% Pd = Diametral pitch of pinion

dP = NP/Pd;
%%
% V
% d = Operating pitch diameter of pinion
% n = Speed

V = pi*d*n/12;
%%
% Wt = Transmitted Load
% H = Power

Wt = 33000*H/V;
%%
% sigma = bending stress
% Ko = overload factor, see table on pg. 766
% Kv = Dynamic Factor

% B = 0.25*(12-Qv)^(2/3);
% A = 50+56*(1-B)
% Kv = ((A+sqrt(V))/A)^B for V in ft/min
% Kv = ((A+sqrt(200*V))/A)^B for V in m/s

%% Ks Formula (Eq.(a), section 14-30, pg. 759)
% Ks = size factor (if Ks < 1 then use 1)

% Ks = 1 or
% Ks = 1/kb or
% Ks = 1.192*(F*sqrt(Y)/P)^0.0535
%% Km Formula (Eq. 14-30, pg. 759)
% F = Net Face Width of Narrowest Member
% Km = Load-Distribution Factor
% 
% Cpf = F/(10*d) - 0.025  for F <= 1 in
% Cpf = F/(10*d) - 0.0375 + 0.0125*F for 1 < F <= 17 in
% Cpf = F/(10*d) - 0.1109 + 0.0207*F - 0.000338*F^2 for 17 < F <= 40 in
% Cmc = 1 (for uncrowned teeth) or  Cmc = 0.8 (for crowned teeth)
% Cpm = 1 (for straddle-mounted pinion w/ S1/S < 0.175) or = 1.1 (for straddle-mounted pinion w/ S1/S >= 0.175)
% A = ?, B = ?, C = ? see table 14-9 pg. 760
% Cma = A + B*F + C*F^2;
% Ce = 0.8 (for gearing adjusted at assembly, compatibility improved by lapping, or both) or = 1 (for all other conditions)
% Cmf = 1 + Cmc*(Cpf*Cpm+Cma*Ce); 
% Km = Cmf

%% KB Formula (Eq. 14-40, pg. 764
% KB = Rim-Thickness factor
% mB = backup ratio
% tR = rim thickness
% ht = tooth height

% mB = tR/ht;
% KB = 1.6*ln(2.242/mB); for mB < 1.2
% KB = 1 for mB >= 1.2
%% J Formula (Fig. 14-6, pg. 753)
% J = Geometry factor for bending strength
% pN = normal base pitch
% Y = form factor
% Z = length of the line of action in the transverse plane
% mN = pN/(0.95*Z); load-sharing ratio
% J = Y/(Kf*mN);

sigma = Wt*Ko*Kv*Ks*(Pd/F)*(Km*KB/J);

