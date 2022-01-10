% Bismillah
% example 1: relative rientation
clc;
clear;
% interior orientation
xo = 0.008; yo = -0.012; f = 152.14;
% exterior orientation parameters
% first image
omega1 = 1.2; phi1 = 2.3; kappa1 = 5.1;
X01 = 1114; Y01 = 862; Z01 = 1500;
% second image
omega2 = 2.5; phi2 = 2.2; kappa2 = 5.7;
X02 = 1926; Y02 = 904; Z02 = 1490;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tie points
XYZ = [1000, 1000, 200;...
       1420, 980, 210;...
       1790, 1700, 155;...
       1800, 340, 180;...
       1095, 295, 166;...
       930, 1650, 170];
xy1 = zeros(6,2);
xy2 = zeros(6,2);
[xy1(1,:)] = BackProjection(XYZ(1,1), XYZ(1,2), XYZ(1,3), omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);
[xy1(2,:)] = BackProjection(XYZ(2,1), XYZ(2,2), XYZ(2,3), omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);
[xy1(3,:)] = BackProjection(XYZ(3,1), XYZ(3,2), XYZ(3,3), omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);
[xy1(4,:)] = BackProjection(XYZ(4,1), XYZ(4,2), XYZ(4,3), omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);
[xy1(5,:)] = BackProjection(XYZ(5,1), XYZ(5,2), XYZ(5,3), omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);
[xy1(6,:)] = BackProjection(XYZ(6,1), XYZ(6,2), XYZ(6,3), omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);

[xy2(1,:)] = BackProjection(XYZ(1,1), XYZ(1,2), XYZ(1,3), omega2, phi2, kappa2, X02, Y02, Z02, xo, yo, f);
[xy2(2,:)] = BackProjection(XYZ(2,1), XYZ(2,2), XYZ(2,3), omega2, phi2, kappa2, X02, Y02, Z02, xo, yo, f);
[xy2(3,:)] = BackProjection(XYZ(3,1), XYZ(3,2), XYZ(3,3), omega2, phi2, kappa2, X02, Y02, Z02, xo, yo, f);
[xy2(4,:)] = BackProjection(XYZ(4,1), XYZ(4,2), XYZ(4,3), omega2, phi2, kappa2, X02, Y02, Z02, xo, yo, f);
[xy2(5,:)] = BackProjection(XYZ(5,1), XYZ(5,2), XYZ(5,3), omega2, phi2, kappa2, X02, Y02, Z02, xo, yo, f);
[xy2(6,:)] = BackProjection(XYZ(6,1), XYZ(6,2), XYZ(6,3), omega2, phi2, kappa2, X02, Y02, Z02, xo, yo, f);

xy1xy2XYZ_Tie = [xy1, xy2, XYZ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Base = 850;
Height = 1450;
[Orient, XYZ_relative]=Relative_Colinear(xy1, xy2, Base, xo , yo , f );
interior= [f, xo , yo ; f, xo , yo ]/1000;
% [DLT1,DLT2,XYZ,iter,Xcap]=Relative_Colinearity2(xy1/1000, xy2/1000, interior, Base, Height);
% [rel,R,T,iter,SigmaX]=Relative_CP(xy1/1000, xy2/1000, interior, Base);
[Orient_cp]=Relative_Coplanar(xy1, xy2, Base, xo , yo , f );
[Orient_abs, residual] = absolute_orientation_M7(XYZ(2:5,:), XYZ_relative(2:5,:) );