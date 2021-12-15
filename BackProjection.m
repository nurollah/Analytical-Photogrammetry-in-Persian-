% Bismillah 
function [xy]=BackProjection(XA, YA, ZA, omega,phi,kappa, X0, Y0, Z0,xo, yo, f)

M=Rottion_Matrix(omega, phi, kappa, 2);

%%%%%%%%%%%%%%%%%%

xa = xo -f * (M(1,1)*(XA-X0) + M(1,2)*(YA-Y0) + M(1,3)*(ZA-Z0))/(M(3,1)*(XA-X0) + M(3,2)*(YA-Y0) + M(3,3)*(ZA-Z0));
ya = yo -f * (M(2,1)*(XA-X0) + M(2,2)*(YA-Y0) + M(2,3)*(ZA-Z0))/(M(3,1)*(XA-X0) + M(3,2)*(YA-Y0) + M(3,3)*(ZA-Z0));

xy = [xa, ya];