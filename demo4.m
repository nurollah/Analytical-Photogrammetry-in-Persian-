%Bismillah
% example 8: 3D Helmert
clc;
clear;

omega = 10;
phi = 15;
kappa = 50;
X0=-200;
Y0=136;
Z0=189;
M=Rottion_Matrix(omega,phi,kappa,2);
R=M';
uvw=[1354;547;135];
xyz=R*uvw+[X0;Y0;Z0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example 9: 3D Affine
clc;
clear;
%
A=[0.6, -0.7, 0.3, -200;
   0.8, 0.60, 0.20, 136;
   -0.1, 0.3, 0.90, 189;
     0  , 0 ,  0 ,  1 ];
%
uvw1=[1354;547;135;1];
xyz1 = A*uvw1;
 
%%%%%%%%%%%%%%%%%%%
% example 10: 3D projective
clc;
clear;
%
A=[0.6, -0.7, 0.3, -200;
   0.8, 0.60, 0.20, 136;
   -0.1, 0.3, 0.90, 189;
   0.01, 0.02,  0 ,  1];
%
uvw1=[1354;547;135;1];
xyz1 = A*uvw1;
xyz1=xyz1/xyz1(4);
%%%%%%%%%%%%%%%%%%%%%%%%%
% example 11: 8-parameters affine
clc;
clear;
A= [0.3,-0.2,0.0004,8060;
    -0.29, 0.18,0.0005, 9600];
uvw1=[1300;650;169;1];
xy = A*uvw1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example 12: co-linearity 
clc;
clear;

% exterior orientation parameters
omega = 2;
phi = 3;
kappa = 10;
X0 = 1114;
Y0 = 862;
Z0 = 1600;
M=Rottion_Matrix(omega,phi,kappa,2);
% interior orientation parameters
xo = 0.008;
yo = -0.12;
f = 152.14;
%
XA = 1300;
YA = 650;
ZA = 169;
%%%%%%%%%%%%%%%%%%

xa = xo -f * (M(1,1)*(XA-X0) + M(1,2)*(YA-Y0) + M(1,3)*(ZA-Z0))/(M(3,1)*(XA-X0) + M(3,2)*(YA-Y0) + M(3,3)*(ZA-Z0));
ya = yo -f * (M(2,1)*(XA-X0) + M(2,2)*(YA-Y0) + M(2,3)*(ZA-Z0))/(M(3,1)*(XA-X0) + M(3,2)*(YA-Y0) + M(3,3)*(ZA-Z0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example 13: co-linearity 
clc;
clear;

% exterior orientation parameters
omega = 1.5;
phi = 2;
kappa = 4;
X0 = 1160;
Y0 = 792;
Z0 = 1500;
M=Rottion_Matrix(omega,phi,kappa,2);
% interior orientation parameters
xo = 0;
yo = 0;
f = 120;
%
XA = 1300;
YA = 650;
ZA = 169;
%%%%%%%%%%%%%%%%%%

xa = xo -f * (M(1,1)*(XA-X0) + M(1,2)*(YA-Y0) + M(1,3)*(ZA-Z0))/(M(3,1)*(XA-X0) + M(3,2)*(YA-Y0) + M(3,3)*(ZA-Z0));
ya = yo -f * (M(2,1)*(XA-X0) + M(2,2)*(YA-Y0) + M(2,3)*(ZA-Z0))/(M(3,1)*(XA-X0) + M(3,2)*(YA-Y0) + M(3,3)*(ZA-Z0));

% file coordinate system
colpp=3840;
rowpp=6912;
pixelsize = 0.012; % milimeters
col = colpp + xa/pixelsize;
row = rowpp - ya/pixelsize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% example 14: rigorous to DLT
clc;
clear;
% exterior orientation parameters
omega = 2;
phi = 3;
kappa = 10;
X0 = 1114;
Y0 = 862;
Z0 = 1600;
M=Rottion_Matrix(omega,phi,kappa,2);
% interior orientation parameters
xo = 0.008;
yo = -0.12;
f = 152.14;

cx = f;
cy = f;
Q = -1/(M(3,1)*(X0) + M(3,2)*(Y0) + M(3,3)*(Z0));

L1 = (xo*M(3,1) - cx*M(1,1))*Q;
L2 = (xo*M(3,2) - cx*M(1,2))*Q;
L3 = (xo*M(3,3) - cx*M(1,3))*Q;
L4 = xo + (M(1,1)*(X0) + M(1,2)*(Y0) + M(1,3)*(Z0))*Q*cx;

L5 = (yo*M(3,1) - cy*M(2,1))*Q;
L6 = (yo*M(3,2) - cy*M(2,2))*Q;
L7 = (yo*M(3,3) - cy*M(2,3))*Q;
L8 = yo + (M(2,1)*(X0) + M(2,2)*(Y0) + M(2,3)*(Z0))*Q*cy;

L9 = M(3,1)*Q;
L10 = M(3,2)*Q;
L11 = M(3,3)*Q;

