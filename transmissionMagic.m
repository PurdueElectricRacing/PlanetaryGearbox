n=0;
GRmax=0;
numPlanets=3;

for Nsun=10:100
    
    for Nbig= 10:100
        
        for Nsmall = 10:100
            
            NC=Nsun+Nbig;
            Nring=NC+Nsmall;
            
            gPS=Nbig/Nsun;
            gRP=Nring/Nsmall;
            %integer=(Nring+Nsun)/numPlanets;
            GR=(Nbig/Nsun*Nring/Nsmall)+1;
            
            if (Nbig+2)<((Nsun+Nbig)*sind(180/numPlanets))
                
                if GR>GRmax
                    GRmax=GR;
                    bestSun=Nsun;
                    bestBig=Nbig;
                    bestSmall=Nsmall;
                    bestRing=Nring;
                end
            end
            n=n+1;
        end
    end
end

fprintf("max gear ratio: %.1f\n",GRmax);
fprintf("\nGear teeth:\nSun: %i\nBig Planet: %i\nSmall Planet: %i\nRing: %i\n\n",bestSun,bestBig,bestSmall,bestRing);




%agma factors of safety chp 13 shigley
