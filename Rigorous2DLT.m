% Bismillah
function [DLT]=Rigorous2DLT(omega,phi,kappa,X0,Y0,Z0,xo,yo,f)
M=Rottion_Matrix(omega,phi,kappa,2);

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

DLT = [L1, L2, L3, L4;...
       L5, L6, L7, L8;...
       L9, L10, L11, 1];
