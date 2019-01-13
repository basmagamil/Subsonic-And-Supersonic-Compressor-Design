function [Pratio, eff, Mrel1, DFr, utr, Cx2, R, phi2, criteria]=supersonic(x)

alpha1=x(1);
phi1=x(2);
psi=x(3);
Cx1=x(4); 
U22U1=x(5); %Ratio U2/U1
zetar=x(6); %rh/rt rotor
soldr=discretize_sold(x(7)); %solidity rotor
h2cr=x(8); %height to chord ratio rotor
zetas=x(9); %rh/rt stator
solds=discretize_sold(x(10)); %solidity stator
h2cs=x(11); %height to chord ratio stator

%consts
cp=1004.5;
Rconst=287;
G=1.4;

%Givens
Tt1=288.16;
Pt1=1.013e5;
mdot=30;

%Facts
r22r1=U22U1;

%Assumptions
Cx2=Cx1;  %change in phase 2
 
%Velocity Triangle
%Vel Tri 1
C1=Cx1./cosd(alpha1);   %Absolute Inlet Velocity
U1=Cx1./phi1;           %Rotor Speed 1
Cth1=C1.*sind(alpha1);  %Absolute Inlet Radial Velocity
Wth1=U1-Cth1;           %Relative Inlet Radial Velocity
W1=sqrt(Cx1.^2+Wth1.^2);%Relative Inlet Velocity
beta1=asind(Wth1./W1);  %Relative Inlet Angle

%Vel Tri 2
U2=U22U1.*U1;           %Rotor Speed 2
phi2=Cx2./U2;           %Flow Coefficient 2
beta2=atand((1-(psi+(1-phi1.*tand(beta1))).*(1./(r22r1.^2)))./phi2);    %Relative Outlet Angle
Wth2=Cx2.*tand(beta2);  %Relative Outlet Radial Velocity
W2=sqrt(Wth2.^2+Cx2.^2);%Relative Outlet Velocity
Cth2=U2-Wth2;           %Absolute Outlet Radial Velocity
C2=sqrt(Cx2.^2+Cth2.^2);%Absolute Outlet Velocity
alpha2=asind(Cth2./C2); %Absolute Outlet Angle

T1=Tt1-(C1.^2)/(2.*cp);     %Inlet Static Temperature
M1=C1./sqrt(G.*Rconst.*T1); %Inlet Mach Number
MFP1=sqrt(G./Rconst).*M1.*(1+((G-1)./2).*(M1.^2)).^(-(G+1)./(2.*(G-1))); %Mass Flow Parameter 1
A1=(mdot.*sqrt(Tt1))./(MFP1.*Pt1.*cosd(alpha1));  %Area 1

%Incidence & Carter relations
%loading graphs
%i
load ('io4.mat'); beta1_io4=io4(:,1); io10_io4=io4(:,2); f_io4=fit(beta1_io4,io10_io4,'smoothingspline');
load ('io5.mat'); beta1_io5=io5(:,1); io10_io5=io5(:,2); f_io5=fit(beta1_io5,io10_io5,'smoothingspline');
load ('io6.mat'); beta1_io6=io6(:,1); io10_io6=io6(:,2); f_io6=fit(beta1_io6,io10_io6,'smoothingspline');
load ('io8.mat'); beta1_io8=io8(:,1); io10_io8=io8(:,2); f_io8=fit(beta1_io8,io10_io8,'smoothingspline');
load ('io10.mat'); beta1_io10=io10(:,1); io10_io10=io10(:,2); f_io10=fit(beta1_io10,io10_io10,'smoothingspline');
load ('io12.mat'); beta1_io12=io12(:,1); io10_io12=io12(:,2); f_io12=fit(beta1_io12,io10_io12,'smoothingspline');
load ('io14.mat'); beta1_io14=io14(:,1); io10_io14=io14(:,2); f_io14=fit(beta1_io14,io10_io14,'smoothingspline');
load ('io15.mat'); beta1_io15=io15(:,1); io10_io15=io15(:,2); f_io15=fit(beta1_io15,io10_io15,'smoothingspline');
load ('io16.mat'); beta1_io16=io16(:,1); io10_io16=io16(:,2); f_io16=fit(beta1_io16,io10_io16,'smoothingspline');
load ('io18.mat'); beta1_io18=io18(:,1); io10_io18=io18(:,2); f_io18=fit(beta1_io18,io10_io18,'smoothingspline');
load ('io20.mat'); beta1_io20=io20(:,1); io10_io20=io20(:,2); f_io20=fit(beta1_io20,io10_io20,'smoothingspline');
%n
load ('n4.mat'); beta1_n4=n4(:,1); n_io4=n4(:,2); f_n4=fit(beta1_n4,n_io4,'smoothingspline');
load ('n5.mat'); beta1_n5=n5(:,1); n_io5=n5(:,2); f_n5=fit(beta1_n5,n_io5,'smoothingspline');
load ('n6.mat'); beta1_n6=n6(:,1); n_io6=n6(:,2); f_n6=fit(beta1_n6,n_io6,'smoothingspline');
load ('n8.mat'); beta1_n8=n8(:,1); n_io8=n8(:,2); f_n8=fit(beta1_n8,n_io8,'smoothingspline');
load ('n10.mat'); beta1_n10=n10(:,1); n_io10=n10(:,2); f_n10=fit(beta1_n10,n_io10,'smoothingspline');
load ('n12.mat'); beta1_n12=n12(:,1); n_io12=n12(:,2); f_n12=fit(beta1_n12,n_io12,'smoothingspline');
load ('n14.mat'); beta1_n14=n14(:,1); n_io14=n14(:,2); f_n14=fit(beta1_n14,n_io14,'smoothingspline');
load ('n15.mat'); beta1_n15=n15(:,1); n_io15=n15(:,2); f_n15=fit(beta1_n15,n_io15,'smoothingspline');
load ('n16.mat'); beta1_n16=n16(:,1); n_io16=n16(:,2); f_n16=fit(beta1_n16,n_io16,'smoothingspline');
load ('n18.mat'); beta1_n18=n18(:,1); n_io18=n18(:,2); f_n18=fit(beta1_n18,n_io18,'smoothingspline');
load ('n20.mat'); beta1_n20=n20(:,1); n_io20=n20(:,2); f_n20=fit(beta1_n20,n_io20,'smoothingspline');
error=0.01;
if (soldr-0.4)<error
    i010=f_io4(beta1); n=f_n4(beta1);
elseif (soldr-0.5)<error
    i010=f_io5(beta1); n=f_n5(beta1); 
elseif (soldr-0.6)<error
    i010=f_io6(beta1); n=f_n6(beta1);
elseif (soldr-0.8)<error
    i010=f_io8(beta1); n=f_n8(beta1);
elseif (soldr-1)<error
    i010=f_io10(beta1); n=f_n10(beta1);
elseif (soldr-1.2)<error
    i010=f_io12(beta1); n=f_n12(beta1);  
elseif (soldr-1.4)<error
    i010=f_io14(beta1); n=f_n14(beta1);
elseif (soldr-1.5)<error
    i010=f_io15(beta1); n=f_n15(beta1);
elseif (soldr-1.6)<error
    i010=f_io16(beta1); n=f_n16(beta1);
elseif (soldr-1.8)<error
    i010=f_io18(beta1); n=f_n18(beta1);
elseif (soldr-2.0)<error
    i010=f_io20(beta1); n=f_n20(beta1);
end

kish=0.7;
kit=0.7;
i0=i010*kish*kit;

a_c=0.5; %circular arc
m=0.23.*((2.*a_c).^2)+beta2./500;
kk=m./((soldr).^0.5);
A=[(n+1),-n;kk,(1-kk)];
b=[beta1-i0;beta2];
XX=A\b;
beta1d=XX(1);
beta2d=XX(2);
i_r=beta1-beta1d;
dev_r=beta2-beta2d;

%Diffusion factor & losses
%DFr=1-(W2/W1)+(abs(r2*Wth2-r1*Wth1))/((2*sold*((r1+r2)/2)*W1));
DFr=1-(W2./W1)+(abs(r22r1.*Wth2-Wth1))./((2.*soldr.*((r22r1./2)+1).*W1));
load ('RotorLoss.mat'); DFr_RL=RotorLoss(:,1); wr3d_RL=RotorLoss(:,2); f_RL=fit(DFr_RL,wr3d_RL,'smoothingspline');
load ('RotorLoss10.mat'); DFr_RL10=RotorLoss10(:,1); wr3d_RL10=RotorLoss10(:,2); f_RL10=fit(DFr_RL10,wr3d_RL10,'smoothingspline');
wr3d_sub=0.7.*f_RL(DFr)+0.3.*f_RL10(DFr);
wr2d_sub=wr3d_sub.*((2.*soldr)./(cosd(beta2d)));

%Transition
P1=Pt1./((Tt1./T1).^(G./(G-1)));
Ttrel1=T1+(W1.^2)./(2.*cp);
Mrel1=W1./sqrt(G.*Rconst.*T1);
Ttrel2=Ttrel1+(U2.^2-U1.^2)./(2.*cp);
Ptrel1=P1.*(Ttrel1./T1).^(G./(G-1));
Ptrel2d=Ptrel1.*(Ttrel2./Ttrel1).^(G./(G-1));

Ptrel2_sub=Ptrel2d-wr2d_sub.*(Ptrel1-P1); 
pi_shock=NSW_totalP(Mrel1);
Ptrel2_super=pi_shock.*Ptrel2_sub;

T2 = Ttrel2 - (W2.^2)./(2.*cp);
Mrel2 = W2./sqrt(G.*Rconst.*T2);
MFPrel2=sqrt(G./Rconst).*Mrel2.*(1+((G-1)./2).*(Mrel2.^2)).^(-(G+1)./(2.*(G-1)));
A2=(mdot.*sqrt(Ttrel2))./(MFPrel2.*Ptrel2_super.*cosd(beta2));

Amr = (A1+A2)./2 ;      %Area Mean
rmr=sqrt((Amr.*(1+zetar))./(4.*pi.*(1-zetar)));%Radius Mean
r2r=2.*rmr./(1+1./r22r1);%Radiua Rotor
r1=r2r./r22r1;           %Radius 1

P2=Ptrel2_super./((1+((G-1)./2).*Mrel2.^2)).^(G./(G-1));
M2=C2./sqrt(G.*Rconst.*T2);
Tt2=T2.*(1+((G-1)./2).*M2.^2);
Pt2=P2.*(1+((G-1)./2).*M2.^2).^(G./(G-1));

%%Stator
%Facts
Tt3=Tt2;
%Assumption: repeating stage
alpha3=alpha1; 
C3=C1;
Cth3=Cth1;

if (solds-0.4)<error   
    i010s=f_io4(alpha2); ns=f_n4(alpha2);
elseif (solds-0.5)<error
    i010s=f_io5(alpha2); ns=f_n5(alpha2); 
elseif (solds-0.6)<error
    i010s=f_io6(alpha2); ns=f_n6(alpha2);
elseif (solds-0.8)<error
    i010s=f_io8(alpha2); ns=f_n8(alpha2);
elseif (solds-1)<error
    i010s=f_io10(alpha2); ns=f_n10(alpha2);
elseif (solds-1.2)<error
    i010s=f_io12(alpha2); ns=f_n12(alpha2);  
elseif (solds-1.4)<error
    i010s=f_io14(alpha2); ns=f_n14(alpha2);
elseif (solds-1.5)<error
    i010s=f_io15(alpha2); ns=f_n15(alpha2);
elseif (solds-1.6)<error
    i010s=f_io16(alpha2); ns=f_n16(alpha2);
elseif (solds-1.8)<error
    i010s=f_io18(alpha2); ns=f_n18(alpha2);
elseif (solds-2.0)<error
    i010s=f_io20(alpha2); ns=f_n20(alpha2);
end

i0s=i010s*kish*kit;

a_c=0.5; %circular arc
m=0.23.*((2.*a_c).^2)+alpha3./500;
kk=m./((solds).^0.5);
A=[(ns+1),-ns;kk,(1-kk)];
b=[alpha2-i0s;alpha3];
XX=A\b;
alpha2d=XX(1);
alpha3d=XX(2);
i_s=alpha2-alpha2d;
dev_s=alpha3-alpha3d;
alpha1d=alpha3d;

staggerR=(beta1d+beta2d)./2;
staggerS=(alpha2d+alpha3d)./2;

%solving for diffusion factor stator
DFsassumed=0.4:0.001:0.47;
tol=0.001;
load ('RotorLoss.mat'); DFs_SL=RotorLoss(:,1); wscos_SL=RotorLoss(:,2); f_SL=fit(DFs_SL,wscos_SL,'smoothingspline');
for j=1:length(DFsassumed)
    wscos=f_SL(DFsassumed(j));
    ws=wscos.*((2.*solds)./(cosd(alpha3d)));  %Stator Loss
    Pt3=Pt2-ws.*(Pt2-P2);                     %Total Pressure 3
    T3 = Tt3 - (C3.^2 /(2.*cp));              %Static Temperature 3
    M3 = C3./sqrt(G.*Rconst.*T3);             %Mach Number 3
    MFP3=sqrt(G./Rconst).*M3.*(1+((G-1)./2).*(M3.^2)).^(-(G+1)./(2.*(G-1)));%Mass Flow Parameter
    A3=(mdot.*sqrt(Tt3))./(MFP3.*Pt3.*cosd(alpha3));  %Area 3
    Ams =(A2+A3)./2;                                  %Area
    rms=sqrt((Ams.*(1+zetas))./(4.*pi.*(1-zetas)));   %Radius
    r2s=r2r+(rms-r2r)/(2.*cosd(staggerS)+1);          %Radius Stator 2
    ratio=2.*(cosd(staggerS)+0.25)/(cosd(staggerS)+0.5);
    r3s=r2r+ratio.*(rms-r2r);                         %Radius Stator 3
    DFscheck=1-C3./C2+(abs(r3s.*Cth3-r2s.*Cth2))./(2.*((r2s+r3s)./2).*solds.*C2);
    %Check 
    if abs(DFsassumed(j)-DFscheck) < tol
        break;
    end
end

P3=Pt3./((1+((G-1)./2).*M3.^2)).^(G./(G-1));

%Calculating pressure ratio and efficiency
Pratio=Pt3./Pt1;    %Stage Pressure Ratio
Tratio=Tt3./Tt1;    %Stage Temperature Ratio
eff=((Pratio.^((G-1)./G))-1)./(Tratio-1);

%calculating no. of blades for 
%rotor
rtr = 2.*rmr./(1+zetar); %rtip rotor
rhr = 2.*rmr-rtr; %rhub rotor
chr=(rtr-rhr)./h2cr; %chord rotor
Sr=chr./soldr; %spacing rotor
Nr=(2.*pi.*rmr)./Sr; %Number of blades in rotor
%stator
rts = 2.*rms./(1+zetas); %rtip stator
rhs = 2.*rms-rts; %rhub stator
chs=(rts-rhs)./h2cs; %chord stator
Ss=chs./solds; %spacing stator
Ns=(2.*pi.*rms)./Ss; %Number of blades in stator

%calculating tip velocities for rotor
omega = U1./r1;    %Radial Velocity
utr = omega.* rtr; %Tip Velocity

R=(((W1.^2)./2)-((W2.^2)./2)+(((U2.^2)-(U1.^2))./2))./((U2.^2)-(U1.^2)+(U1.*Wth1)-(U2.*Wth2));

criteria=min(-(0.6.*(Pratio)+0.4.*(eff)));
if((Mrel1>1)&&(Mrel1<=1.25)&&(DFr>=0.3)&&(DFr<=0.47)&&(utr<=400)&&(Cx2>=100)&&(Cx2<=150)&&(R<=1)&&(R>=0.5)&&(phi2>=0.25)&&(phi2<=0.5)) 
 % Satisfied 
    criteria=min(-(0.6.*(Pratio)+0.4.*(eff)));
else
    criteria=455567;
end
end
