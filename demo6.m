% Bismillah
% example 1: 8-parameters affine
clc;
clear;
%format long;
A= [0.29,-0.21,0.0004,8860;
    -0.30, 0.18,0.0005, 9600];

XYZ1=[1300, 900,  1000 , 1800, 800;...
      650,  1250, 1860,  890,  1600;...
      169,  120,  210,   245,  100;...
      1,    1,    1,     1,     1];
xy = A*XYZ1;
%
xy = xy';
xy = round(xy);
XYZ1 = XYZ1';

[Par, res]=Affine_3D_2D(XYZ1, xy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example 2: space resection using DLT
clc;
clear;
% interior orientation
xo = 0.008; yo = -0.012; f = 152.14;
% exterior orientation parameters
% first image
% omega1 = 3.1; phi1 = 2.3; kappa1 = 4.2;
% X01 = 1114; Y01 = 862; Z01 = 1600;
omega1 = 1.1; phi1 = 3.3; kappa1 = 6.2;
X01 = 1120; Y01 = 962; Z01 = 1550;
XYZ =[1300, 900,  1000 , 1800, 800, 1110;...
      650,  1250, 1860,  890,  1600, 987;...
      169,  120,  210,   245,  100,  251];
 XYZ = XYZ';
 xy1 = zeros(6,2); 
[xy1(1,:)] = BackProjection(XYZ(1,1), XYZ(1,2), XYZ(1,3), omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);
[xy1(2,:)] = BackProjection(XYZ(2,1), XYZ(2,2), XYZ(2,3), omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);
[xy1(3,:)] = BackProjection(XYZ(3,1), XYZ(3,2), XYZ(3,3), omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);
[xy1(4,:)] = BackProjection(XYZ(4,1), XYZ(4,2), XYZ(4,3), omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);
[xy1(5,:)] = BackProjection(XYZ(5,1), XYZ(5,2), XYZ(5,3), omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);
[xy1(6,:)] = BackProjection(XYZ(6,1), XYZ(6,2), XYZ(6,3), omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);
%
[DLT1]=Rigorous2DLT(omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);

xy = round(xy1,4);
[Par_DLT, residuals]=DLT_Resection(XYZ, xy);
Par_DLT(12,1)=1;
DLT = reshape(Par_DLT,4,3);
DLT = DLT';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
% interior orientation
xo = 0.008; yo = -0.012; f = 152.14;
% exterior orientation parameters
% first image
omega1 = 2; phi1 = 3; kappa1 = 6.1;
X01 = 1114; Y01 = 862; Z01 = 1600;
% coordinates of an optional ground point 
XYZ = [1260, 1410, 210;...
       1000, 1650, 150;...
       1850, 900, 100;...
       900, 1180, 180];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process
xy1 = zeros(4,2);
[xy1(1,:)] = BackProjection(XYZ(1,1), XYZ(1,2), XYZ(1,3), omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);
[xy1(2,:)] = BackProjection(XYZ(2,1), XYZ(2,2), XYZ(2,3), omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);
[xy1(3,:)] = BackProjection(XYZ(3,1), XYZ(3,2), XYZ(3,3), omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);
[xy1(4,:)] = BackProjection(XYZ(4,1), XYZ(4,2), XYZ(4,3), omega1, phi1, kappa1, X01, Y01, Z01, xo, yo, f);

xy1=round(xy1,4);
% xy1 = xy1 + rand(4,2)/100-0.005;% manually adding residuals  (these residuals are less than 10 micron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Exterior]=Space_Resection(XYZ, xy1, xo , yo , f );

