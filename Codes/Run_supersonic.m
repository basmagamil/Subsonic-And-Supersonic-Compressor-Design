clc;clear all;close all;
Nv = 11;
tic 
%xbound=[alpha1, phi1, psi, Cx1, U22U1, zetar, soldr, h2cr, zetas, solds, h2cs]
xbound =[0 15; %alpha1
    0.25 0.5; %phi1
    0.1 0.4; %psi 
    100 150; %Cx1
    1 1.05; %U22U1
    0.4 0.9; %zetar
    0.4 2; %soldr
    3 5; %h2cr
    0.6 0.8; %zetas
    0.4 2; %solds
    3 5]; %h2cs 
X0 = xbound(:,1) + (xbound(:,2) - xbound(:,1)).*rand(Nv,1);
Lb = xbound(:,1);
Ub = xbound(:,2);
options = gaoptimset('TimeLimit',12000,'PopulationSize',300,'Generations',800, 'TolCon', 1e-5,...
'TolFun', 1e-5,'PlotFcns',@gaplotbestf);
options.Display='iter';
[x, fval, exitflag,output] = ga(@optimize_supersonic,11, [], [], [], [], Lb, Ub, [],options);
[Pratio, eff, Mrel1, DFr, utr, Cx2, R, phi2, criteria]=supersonic(x)
x=[x(1) x(2) x(3) x(4) x(5) x(6) discretize_sold(x(7)) x(8) x(9) discretize_sold(x(10)) x(11)]
criteria= fval
toc
