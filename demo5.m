% Bismillah
% this demo is used for fivth chapter of analytical photogrammetry
% this chapter focouses on space intersection
% Nurollah Tatar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example 1: space intersection using DLT
clc;
clear;
% interior orientation
xo = 0.008; yo = -0.012; f = 152.14;
% exterior orientation parameters
% first image
omega1 = 2; phi1 = 3; kappa1 = 6.1;
X01 = 1114; Y01 = 862; Z01 = 1600;
% second image
omega2 = 3.5; phi2 = 2; kappa2 = 4.7;
X02 = 1966; Y02 = 904; Z02 = 1590;

% coordinates of an optional ground point 
XA = 1260; 
YA = 1410; 
ZA = 210;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process
[xy1]=BackProjection(XA, YA, ZA, omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);
[xy2]=BackProjection(XA, YA, ZA, omega2, phi2, kappa2, X02, Y02, Z02, xo, yo, f);

[DLT1]=Rigorous2DLT(omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);
[DLT2]=Rigorous2DLT(omega2, phi2, kappa2, X02, Y02, Z02, xo, yo, f);
%%%%%%%%%%%%%%%%%%
[XYZ]=Intersection_By_DLT(xy1, xy2, DLT1, DLT2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example 2: space intersection using co linear equations
exterior = [omega1, phi1, kappa1, X01, Y01, Z01;...
            omega2, phi2, kappa2, X02, Y02, Z02];
interior = [xo, yo, f];
[XYZ2]=Intersection_by_CoLinear(xy1, xy2, exterior, interior);






