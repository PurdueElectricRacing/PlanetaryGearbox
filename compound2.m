function compound2(gr,planets, minOD, maxOD, sunOD_min)
%% NOTE
%This function currently produces potential coumpound planetary gear
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
%minOD = 2 * (sun radius + big planet diameter) (minimum value)
%maxOD = 2 * (sun radius + big planet diameter) (maximum value)
%sunOD_min = minimum OD of sun gear that allows a motor shaft to properly
%fit through

pd_stand = [16 18 20 22 24]; %diametral pitch

%% CALCULATIONS
for ring_OD = minOD:0.01:maxOD %ring_OD is imaginary ring gear around first stage
    
  %%first stage calculations
    Nr = zeros(1, numel(pd_stand));
    %Finds all the possible ring gear teeth counts for specific desired OD
    %and diametral pitch value
    for i=1:numel(pd_stand)
        Nr(i) = round(pd_stand(i)*ring_OD);
    end
    
    %Finds sun and planet gears
    Ns = round(pd_stand.*sunOD_min);
    Np = zeros(1,numel(Ns));
    for i=1:numel(Nr)
        real = (Nr(i)+Ns(i))/planets; %condition 2
        while floor(real) ~= real
            Ns(i) = Ns(i) + 1;
            real = (Nr(i)+Ns(i))/planets;
        end
        planet = (Nr(i)-Ns(i))/2; %condition 1
        if floor(planet) == planet
            Np(i) = (Nr(i)-Ns(i))./2;
        else
            Np(i) = 0;
        end
    end
    
  %%second stage calculations
    ring_OD2 = 0.8*ring_OD;
    Nr2 = round(ring_OD2*pd_stand);
    Np2 = zeros(1,numel(Nr2));
    Ns2 = zeros(1,numel(Nr2));
    for i=1:numel(Np)
       if Np(i) ~= 0
           Np2(i) = round(1 ./ ((gr - 1) .* Ns(i) ./ Np(i) ./ Nr2(i))); 
           Ns2(i) = Nr2(i) - 2*Np2(i); %condition 1
           while floor(Ns2) ~= Ns2
               Np2(i) = Np2(i) + 1;
               Ns2(i) = Nr2(i) - 2*Np2(i);
           end
       else
           Np2(i) = 0;
       end
    end
    
    fprintf('Ring OD: %.2f', ring_OD)
    pd_2 = transpose(pd_stand);
    Ns_t = transpose(Ns);
    Np_t = transpose(Np);
    Np2_t = transpose(Np2);
    Nr2_t = transpose(Nr2);
    Final = cat(2,pd_2, Ns_t,Np_t,Np2_t,Nr2_t);
    colNames = {'P_d','Ns','Np_big','Np_small','Nr'};
    gear_table = array2table(Final, 'VariableNames',colNames)
    
end