function GearCalculator(gr,range,planets, ringOD_min, ringOD_max)
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
pd_stand = [14 16 18 20 22]; %diametral pitch

%% CALCULATIONS
for ring_OD = ringOD_min:0.1:ringOD_max
    Nr = zeros(1, numel(pd_stand));
    %Finds all the possible ring gear teeth counts for specific desired OD
    %and diametral pitch value
    for i=1:numel(pd_stand)
        Nr(i) = round(pd_stand(i)*ring_OD);
    end
    
    %Finds all the sun gear teeth counts that allow to have desired gear
    %ratio for each diametral pitch value
    Ns = zeros(5,2);
    for i=1:numel(Nr)
        count = 1;
        for j=1:25
            gr_init = 1 + Nr(i)./j;
            if gr_init < (gr-range)
                break
            elseif gr_init >= (gr - range) && gr_init <= (gr+range)
                Ns(i,count) = j;
                count = count + 1;
            end  
        end
    end

    %Finds all possible planet teeth counts that can actually mesh with
    %previously found sun and ring gears
    gearR = zeros(size(Ns));
    Np = zeros(size(Ns));
    [row,col] = size(Ns);
    for i=1:numel(Nr)
        for j=1:col
            if Ns(i,j) ~= 0
                real = (Nr(i)+Ns(i,j))/planets;
                if floor(real) == real
                    planet = (Nr(i)-Ns(i,j))/2;
                    if floor(planet) == planet
                        Np(i,j) = (Nr(i)-Ns(i,j))./2;
                        gearR(i,j) = 1 + Nr(i)./Ns(i,j);
                    else
                        Np(i,j) = 0;
                        gearR(i,j) = 0;
                    end
                else 
                    Np(i,j) = 0;
                    gearR(i,j) = 0;
                end
            else
                break
            end
        end
    end
    
    %output all the possible setups for a set ring gear OD
    fprintf('Ring OD: %.2f', ring_OD)
    Nr_2 = transpose(Nr);
    pd_2 = transpose(pd_stand);
    Final = cat(2,pd_2,Nr_2,Ns(:,1),Np(:,1),Ns(:,2),Np(:,2),gearR(:,1), gearR(:,2));
    colNames = {'P_d','Nr','Ns','Np','Ns2','Np2','GR1','GR2'};
    gear_table = array2table(Final, 'VariableNames',colNames)
    
    
end

    