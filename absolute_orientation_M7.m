% Bismillah
function [Orient, residual]=absolute_orientation_M7(XYZ, xyz )
% about function: this function is used to do compute parameters of absolute orientation 
% based on 3D Helmert method.
% this matlab code implemented by Nurollah Tatar, PhD  in photogrammetry at
% University of Theran, Tehran, Iran. Email: n.tatar@ut.ac.ir
% inputs:
% XYZ= [X,Y,Z] the coordinate in object space (meters). this points must be more than 3 pts
% xyz= [x,y,z] the coordinate in model space (meters). this points must be more than 3 pts
% compute initial values 
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
n = size(xyz,1);
%%%%%%%%%%%%%%%%%%%%%%%
% compute initial values for exterior orientation
A = zeros(2*n,4);
A(1:2:end,1) = x;
A(1:2:end,2) = y;
A(2:2:end,1) = y;
A(2:2:end,2) = -x;
A(1:2:end,3) = 1;
A(2:2:end,4) = 1;
%
L = zeros(2*n,1);
L(1:2:end,1) = XYZ(:,1);
L(2:2:end,1) = XYZ(:,2);
XX = inv(A'*A)*A'*L;
a = XX(1);
b = XX(2);
c = XX(3);
d = XX(4);
landa0 = sqrt(a^2 + b^2);
X0 = c;
Y0 = d;
Z0 =  mean(XYZ(:,3)) - landa0*mean(z);
kapa0 = atan(-b/a);
omega0 = 0;
phi0 = 0;
%
for i=1:40
    % the angles are according to degrees
    Ck = cos(kapa0);
    Sk = sin(kapa0);
    Co = cos(omega0);
    So = sin(omega0);
    Cph = cos(phi0);
    Sph = sin(phi0);
    
     Mx = [1, 0, 0; 0, Co, So; 0, -So, Co];
     My = [Cph, 0, -Sph; 0, 1, 0; Sph, 0, Cph];
     Mz = [Ck, Sk, 0; -Sk, Ck, 0; 0, 0, 1];
     M = Mz*My*Mx;
     R = M';
     % 
     r = R(1,1)*x + R(1,2)*y + R(1,3)*z;
     s = R(2,1)*x + R(2,2)*y + R(2,3)*z;
     q = R(3,1)*x + R(3,2)*y + R(3,3)*z;
     
     roundF_omega = zeros(n,1);
     roundF_phi = landa0*(-Sph*Ck*x + Sph*Sk*y + Cph*z);
     roundF_kapa = landa0*(R(1,2)*x - R(1,1)*y);
     roundF_landa = r;
     roundF_X0 = ones(n,1);
     %
     roundG_omega =  -landa0*q;
     roundG_phi =  landa0*So*r;
     roundG_kapa =  landa0*(R(2,2)*x - R(2,1)*y);
     roundG_landa =  s;
     roundG_Y0 =  ones(n,1);
     % 
     roundH_omega =  s*landa0;
     roundH_phi =  -r*landa0*Co;
     roundH_kapa =  landa0*(R(3,2)*x - R(3,1)*y);
     roundH_landa =  q;
     roundH_Z0 =  ones(n,1);     
     %     
     FX0 = landa0*r + X0;
     GX0 = landa0*s + Y0;
     HX0 = landa0*q + Z0;
     %
     A = zeros(3*n,7);
     A(1:3:end,1:5) = [roundF_omega, roundF_phi, roundF_kapa, roundF_landa, roundF_X0];
     A(2:3:end,1:4) = [roundG_omega, roundG_phi, roundG_kapa, roundG_landa];
     A(2:3:end,6) = roundG_Y0;
     A(3:3:end,1:4) = [roundH_omega, roundH_phi, roundH_kapa, roundH_landa];
     A(3:3:end,7) = roundH_Z0;
     %
     L = zeros(3*n,1);
     L(1:3:end,1) =  XYZ(:,1) - FX0;
     L(2:3:end,1) =  XYZ(:,2) - GX0;
     L(3:3:end,1) =  XYZ(:,3) - HX0;
     
     XX = inv(A'*A)*A'*L;%
     dx1 = max(abs(XX));
   
   if dx1 <0.0000001
       break;
   else
       omega0 = omega0 + XX(1);
       phi0 = phi0 + XX(2);
       kapa0 = kapa0 + XX(3); 
       landa0 = landa0 + XX(4);
       X0 = X0 + XX(5);
       Y0 = Y0 + XX(6);
       Z0 = Z0 + XX(7);
   end             
     
end
residual = XYZ'-landa0*R*xyz' - repmat([X0;Y0;Z0],1,n);
residual = residual';
% convert radian to degree
omega = omega0 * 180/pi();
phi = phi0 * 180/pi();
kapa = kapa0 * 180/pi();

Orient = [omega; phi; kapa; landa0; X0; Y0; Z0];