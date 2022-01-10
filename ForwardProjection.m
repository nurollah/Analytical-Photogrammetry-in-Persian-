function [XY]=ForwardProjection(xa, ya, ZA, omega,phi,kappa, X0, Y0, Z0,xo, yo, f)

M=Rottion_Matrix(omega, phi, kappa, 2);
R = M';
%%%%%%%%%%%%%%%%%%

XA = X0 +(ZA-Z0) * (R(1,1)*(xa-xo) + R(1,2)*(ya-yo) + R(1,3)*(-f))./(R(3,1)*(xa-xo) + R(3,2)*(ya-yo) + R(3,3)*(-f));
YA = Y0 +(ZA-Z0) * (R(2,1)*(xa-xo) + R(2,2)*(ya-yo) + R(2,3)*(-f))./(R(3,1)*(xa-xo) + R(3,2)*(ya-yo) + R(3,3)*(-f));

XY = [XA, YA];