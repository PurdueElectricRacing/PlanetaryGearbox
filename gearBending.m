
%% Geometry
dg=5; %placeholder 
n=1800; %rpm
H=4; %hp placeholder
Y=.966; %Lewis factor changes w/ #teeth
Np=17; %guess
Pd=10; %transverse diametral pitch
dp=Np/Pd; %pinion
F=1.5; %face width (gear thickness)
Qv=6; %3-7, depends on gear

%% Equations
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




%% final eq
bending=Wt*Ko*Kv*Ks*Pd/F*Km*Kb/J;


