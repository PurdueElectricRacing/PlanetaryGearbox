function GearCalculator(gr,range,planets, sunOD_min, maxOD)
%% NOTE
%This function currently produces potential single stage planetary gear
%setups for different diametral pitch values (defined below) and other
%input parameters. 
%
%There is currently no excel sheet output. Rather the output is in the
%following format:
%%%
%ring OD:
%Nr Ns Np Ns2 Np2 (for pd = lowest value)
%.  .  .  .   .
%Nr Ns Np Ns2 Np2 (for pd = highest value)
%
%Where N = number of teeth. Note asingle ring gear setup may have multiple 
%planet and sun setups, hence Ns2 and Np2)
%% INPUTS
%gr = desired gear ratio
pd_stand = [16 18 20 22 24]; %diametral pitch
motorS_d = .709;
pA = 20; %pressure angle

%% CALCULATIONS
Ns = round(sunOD_min.*pd_stand);
sOD_new = Ns./pd_stand;

pB_OD = maxOD - sOD_new/2;
NpB = round(pB_OD*pd_stand);

for i = 1:numel(NpB)
    real = NpB(i)

ring_OD = 0.8*maxOD;
ring_teeth = ring_OD*pd_stand;

    
end

    

    