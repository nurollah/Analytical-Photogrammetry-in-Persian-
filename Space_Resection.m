% Bismillah
function [Exterior]=Space_Resection(XYZ, xy, xo , yo , f )
% about function: this function is used to do compute exterior orientation from
% ground control points based on space resection by co-linear equation.
% this matlab code implemented by Nurollah Tatar, PhD  in photogrammetry at
% University of Theran, Tehran, Iran. Email: n.tatar@ut.ac.ir
% inputs:
% xy= image coordinates. this points must be more than 3 pts
% XYZ= [X,Y,Z] the coordinates of ground control points (meters). this points must be more than 3 pts
% xo,yo : coordinate of principle point 
%  f is focal length
%
m =size(XYZ,1);
% compute initial values for exterior orientation

A = zeros(2*m,4);
A(1:2:end,1) = xy(:,1)/1000;% milimeters to meters
A(1:2:end,2) = xy(:,2)/1000;% milimeters to meters
A(2:2:end,1) = xy(:,2)/1000;% milimeters to meters
A(2:2:end,2) = -xy(:,1)/1000;% milimeters to meters
A(1:2:end,3) = 1;
A(2:2:end,4) = 1;
%
L = zeros(2*m,1);
L(1:2:end,1) = XYZ(:,1);
L(2:2:end,1) = XYZ(:,2);
XX = inv(A'*A)*A'*L;
a = XX(1);
b = XX(2);
c = XX(3);
d = XX(4);
landa = sqrt(a^2 + b^2);
X0 = c;
Y0 = d;
Z0 = f*landa/1000 + mean(XYZ(:,3));
kapa0 = atan(-b/a);
omega0 = 0;
phi0 = 0;

for i=1:20
    
     deltaX = XYZ(:,1) - X0;
     deltaY = XYZ(:,2) - Y0;
     deltaZ = XYZ(:,3) - Z0;
    % the angles are according to degrees
    Ck = cos(kapa0);
    Sk = sin(kapa0);
    Co = cos(omega0);
    So = sin(omega0);
    Cph = cos(phi0);
    Sph = sin(phi0);
%     M = [Ck*Cph       Ck*Sph*So+Sk*Co   -Ck*Sph*Co+Sk*So;...
%       -Sk*Cph      -Sk*Sph*So+Ck*Co    Sk*Sph*Co+Ck*So;...
%         Sph            -Cph*So             Cph*Co];
     Mx = [1, 0, 0; 0, Co, So; 0, -So, Co];
     My = [Cph, 0, -Sph; 0, 1, 0; Sph, 0, Cph];
     Mz = [Ck, Sk, 0; -Sk, Ck, 0; 0, 0, 1];
     M = Mz*My*Mx;
     % 
     r = M(1,1)*deltaX + M(1,2)*deltaY + M(1,3)*deltaZ;
     s = M(2,1)*deltaX + M(2,2)*deltaY + M(2,3)*deltaZ;
     q = M(3,1)*deltaX + M(3,2)*deltaY + M(3,3)*deltaZ;
     
     roundF_omega = (f*ones(m,1)./(q.^2)).*(r.*(-M(3,3)*deltaY+M(3,2)*deltaZ) - q.*(-M(1,3)*deltaY+M(1,2)*deltaZ));
     roundF_phi = (f*ones(m,1)./(q.^2)).*(r.*(Cph*deltaX + So*Sph*deltaY - Co*Sph*deltaZ)-...
                    q.*(-Sph*Ck*deltaX + So*Cph*Ck*deltaY -Co*Cph*Ck*deltaZ ));
     roundF_kapa = (-f*ones(m,1)./q).*(M(2,1)*deltaX + M(2,2)*deltaY + M(2,3)*deltaZ);
     
     roundF_X0  = (-f*ones(m,1)./(q.^2)).*(r*M(3,1) - q*M(1,1));
     roundF_Y0  = (-f*ones(m,1)./(q.^2)).*(r*M(3,2) - q*M(1,2));
     roundF_Z0  = (-f*ones(m,1)./(q.^2)).*(r*M(3,3) - q*M(1,3));
     %
     roundG_omega = (f*ones(m,1)./(q.^2)).*(s.*(-M(3,3)*deltaY+M(3,2)*deltaZ) - q.*(-M(2,3)*deltaY+M(2,2)*deltaZ));
     roundG_phi = (f*ones(m,1)./(q.^2)).*(s.*(Cph*deltaX + So*Sph*deltaY - Co*Sph*deltaZ)-...
                    q.*(Sph*Sk*deltaX - So*Cph*Sk*deltaY + Co*Cph*Sk*deltaZ ));
     roundG_kapa = (f*ones(m,1)./q).*(M(1,1)*deltaX + M(1,2)*deltaY + M(1,3)*deltaZ);
     roundG_X0  = (-f*ones(m,1)./(q.^2)).*(s*M(3,1) - q*M(2,1));
     roundG_Y0  = (-f*ones(m,1)./(q.^2)).*(s*M(3,2) - q*M(2,2));
     roundG_Z0  = (-f*ones(m,1)./(q.^2)).*(s*M(3,3) - q*M(2,3));
     %
     FX0 = xo -f * r./q;
     GX0 = yo -f * s./q;
     %
     A = zeros(2*m, 6);
     A(1:2:end,1:6)=[roundF_omega, roundF_phi, roundF_kapa, roundF_X0, roundF_Y0, roundF_Z0];
     A(2:2:end,1:6)=[roundG_omega, roundG_phi, roundG_kapa, roundG_X0, roundG_Y0, roundG_Z0];
     %
     L = zeros(2*m, 1);
     L(1:2:end) = xy(:,1)-FX0 + roundF_omega*omega0 + roundF_phi*phi0 +...
                  roundF_kapa*kapa0 + roundF_X0*X0 + roundF_Y0*Y0 + roundF_Z0*Z0;
     L(2:2:end) = xy(:,2)-GX0 + roundG_omega*omega0 + roundG_phi*phi0 +...
                  roundG_kapa*kapa0 + roundG_X0*X0 + roundG_Y0*Y0 + roundG_Z0*Z0;
     %
     XX = inv(A'*A)*A'*L;
   dx = XX-[omega0; phi0; kapa0; X0; Y0; Z0];
   dx = max(abs(dx));
   
   if dx <0.0000001
       break;
   else
       omega0 = XX(1);
       phi0 = XX(2);
       kapa0 = XX(3);
       X0 = XX(4);
       Y0 = XX(5);
       Z0 = XX(6);
   end
              
     
end
omega = XX(1)*180/pi();
phi = XX(2)*180/pi();
kapa = XX(3)*180/pi();
X = XX(4);
Y = XX(5);
Z = XX(6);
Exterior = [omega; phi; kapa; X; Y; Z];