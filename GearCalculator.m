%% INPUTS
% FOS_min
gr = 7.25;
range = 0.25; %+- 0.25 from gr
planet_n = 3;
% planet_OD
% face_w
dp_stand = [6 8 10 12 16]; %diametral pitch
% pa_stand = [20 22.5 25];
% material
% quality

%% CALCULATIONS
for ring_OD = 5:0.1:7
    Nr = zeros(1, numel(dp_stand));
    for i=1:numel(dp_stand)
        Nr(i) = round(dp_stand(i)*ring_OD);
    end

    Ns = zeros(3,2);
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

    Np = zeros(size(Ns));
    [row,col] = size(Ns);
    for i=1:numel(Nr)
        for j=1:col
            if Ns(i,j) ~= 0
                real = (Nr(i)+Ns(i,j))/planet_n;
                if floor(real) == real
                    planet = (Nr(i)-Ns(i,j))/2;
                    if floor(planet) == planet
                        Np(i,j) = (Nr(i)-Ns(i,j))./2;
                    else
                        Np(i,j) = 0;
                    end
                else 
                    Np(i,j) = 0;
                end
            else
                break
            end
        end
    end
    ring_OD
    Nr_2 = transpose(Nr);
    Final = cat(2,Nr_2,Ns(:,1),Np(:,1),Ns(:,2),Np(:,2))
end

    