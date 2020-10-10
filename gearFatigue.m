function [sigC, sigContactAllowed]= gearFatigue(Np,Ng,Pd,d,n,H,Dp,F,Cmc,rg,rp)


%% inputs
Np=1; %teeth pinion
Ng=1; %teeth gear
Pd=1; %pitch diameter
d=1; %gear diameter
n=1; %rpm
H=1; %hp
dp=1; %diameter pinion
F=1; %face widtj (gear thickness)
Cmc=1; %crowned/not
rg=1; %radius gear
rp=1; %radius pinion
%% constants
%assume steel on steel gear [Table 14-8]
Cp=2300;
Ko = 1.25;
Qv=6; %assume commercial quality
Cf=1;
Sc=276000; %table 14-6
Zn=1; %Figure 14-15
CH=1.06; %figure 14-12 depends on mg
Kt=1;
Kr=1; %0.99 reliability
%% equations
V=3.14*d*n/12;
Wt=33000*H/V;

B=.25*(12-Qv)^(2/3);
A=50+56*(1-B);
Kv=((A+sqrt(V)/A)^B);

Ks=1;

Cpf=F/10/dp-.025;
Cpm=1; %maybe 1.1
Cma=.127+.0158*F-(.93*10^(-4))*F^2; %commercial enclosed units
Ce=1; %maybe .8
Km=1+Cmc*(Cpf*Cpm+Cma*Ce);
mg=Ng/Np;
psit=asind(2*(1/dp+1/d)/(1/rg+1/rp));
I=cosd(psit)*sind(psit)/2/(mg-1)*(mg);

sigC=Cp*(Wt*Ko*Kv*Ks*Km*Cf/dp/F/I)^(1/2);
SH=Sc*Zn*CH/(Kt*Kr*sigC);
sigContactAllowed=Sc*Zn*CH/(SH*Kt*Kr);



end

