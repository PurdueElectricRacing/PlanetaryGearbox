n=0;
GRmax=0;
numPlanets=3;
mod=1; %assume
bend=0;
%% Geometry
dg=5; %gear diameter
n=1800; %rpm
T=100; %guess
H=T*n/5252; %hp placeholder
Y=.966; %Lewis factor changes w/ #teeth
mod=1; %assume
F=1; %face width (gear thickness)
Qv=6; %3-7, depends on gear

%%
for Nsun=10:100
      
        for Nplanet = 10:100
            Pd=Nplanet*mod; %transverse diametral pitch
            dp=Nplanet/Pd; %pinion
            Nring=NC+Nplanet;
            
            gPS=Nbig/Nsun;
            gRP=Nring/Nplanet;
            %integer=(Nring+Nsun)/numPlanets;
            GR=(Nbig/Nsun*Nring/Nplanet)+1;
            
            if (Nbig+2)<((Nsun+Nbig)*sind(180/numPlanets))       
                V=3.14*dp*n/12;
                Wt=33000*H/V;
                
                q=.8; %guess
                Kt=2; %guess
                Kf=1+q*(Kt-1);
                
                Ko=1; %assumed uniform moderate shock
                
                B=.25*(12-Qv)^(2/3);
                A=50+56*(1-B);
                Kv=((A+sqrt(V))/A)^B;
                kb=.91*dg^-.157; %.11<=d<=2
                Ks=1/kb;
                Cmc=1; %assume uncrowned .8 if crowned
                Cpf=F/10/dp-.025;
                Cpm=1; %assume S1/s<.175
                Cma=.247+.0167*F+(-.765*10^-4); %assume open gearing
                Ce=1; %other conditions
                Km=1+Cmc*(Cpf*Cpm+Cma*Ce);
                Kb=1; %assume thick enough rim
                J=Y/Kf;
                
     
                bending=Wt*Ko*Kv*Ks*Pd/F*Km*Kb/J;
                %%
                if GR>GRmax
                   
                    bend=bending;
                    GRmax=GR;
                    bestSun=Nsun;
                    bestSmall=Nplanet;
                    bestRing=Nring;
                    
                    
                end
            end
            n=n+1;
    
    end
end

fprintf("max gear ratio: %.1f\n",GRmax);
fprintf("\nGear teeth:\nSun: %i\nBig Planet: %i\nSmall Planet: %i\nRing: %i\n\n",bestSun,bestBig,bestSmall,bestRing);
fprintf("max gear bending: %.1f\n",bend);



%agma factors of safety chp 13 shigley
