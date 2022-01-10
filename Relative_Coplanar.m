% Bismillah
function [Orient]=Relative_Coplanar(xy1, xy2, Base, xo , yo , f )
% about function: this function is used to do compute relative orientation parameters
% based on co-planar equation.
% this matlab code implemented by Nurollah Tatar, PhD  in photogrammetry at
% University of Theran, Tehran, Iran. Email: n.tatar@ut.ac.ir
% inputs:
% xy1= [x1,y1] the coordinate in image 1 (left image). this points must be more than 5 pts
% xy2= [x2,y2] the coordinate in image 2 (rigth image). this points must be more than 5 pts
% xo,yo : coordinate of principle point 
%  f is focal length
% Base: it i base line between stereo camera according to meter. 
xy1 = xy1 / 1000;% milimeters to meters
xy2 = xy2 / 1000;% milimeters to meters
xo = xo/1000;% milimeters to meters
yo = yo/1000;% milimeters to meters
f = f/1000;% milimeters to meters
% compute initial values 
x1 = xy1(:,1);
y1 = xy1(:,2);
x2 = xy2(:,1);
y2 = xy2(:,2);
%%%%%%%%%%%%%%%%%%%%%%%
kapa0 = 0;
omega0 = 0;
phi0 = 0;

n = size(x1,1);
Bx =Base;
By = 0;
Bz = 0;

for i=1:40
    
     deltaX1 = x1 - xo;
     deltaY1 = y1 - yo;
     deltaX2 = x2 - xo;
     deltaY2 = y2 - yo;
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
     r = R(1,1)*deltaX2 + R(1,2)*deltaY2 - R(1,3)*f;
     s = R(2,1)*deltaX2 + R(2,2)*deltaY2 - R(2,3)*f;
     q = R(3,1)*deltaX2 + R(3,2)*deltaY2 - R(3,3)*f;
     
     roundF_By = -f*r - deltaX1.*q;
     roundF_Bz = -deltaY1.*r + deltaX1.*s;
     roundF_omega = (Bx*f + Bz*deltaX1).*(-q) + (Bx*deltaY1 - By*deltaX1).*s;
     roundF_phi = (-By*f - Bz*deltaY1).*(-Sph*Ck*deltaX2 + Sph*Sk*deltaY2 - Cph*f)+...
                   (Bx*f + Bz*deltaX1).*(R(2,2)*deltaX2 - R(2,1)*deltaY2 - Sph*So*f)+...
                   (Bx*deltaY1 - By*deltaX1).*(-Co*r);
     roundF_kapa = (-By*f - Bz*deltaY1).*(R(1,2)*deltaX2 - R(1,1)*deltaY2)+...
                   (Bx*f + Bz*deltaX1).*(R(2,2)*deltaX2 - R(2,1)*deltaY2)+...
                   (Bx*deltaY1 - By*deltaX1).*(R(3,2)*deltaX2 - R(3,1)*deltaY2);
     
     FX0 = (-By*f - Bz*deltaY1).*r +...
           (Bx*f + Bz*deltaX1).*s +...
           (Bx*deltaY1 - By*deltaX1).*q;
     %
     A =[roundF_By, roundF_Bz, roundF_omega,  roundF_phi, roundF_kapa];
     %
     L = -FX0;    
     
     XX = inv(A'*A)*A'*L;%
   dx1 = max(abs(XX));
   
   if dx1 <0.000001
       break;
   else
       By = By + XX(1);
       Bz = Bz + XX(2);
       omega0 = omega0 + XX(3);
       phi0 = phi0 + XX(4);
       kapa0 = kapa0 + XX(5);       
   end             
     
end
By = By + XX(1);
Bz = Bz + XX(2);
omega = omega0 + XX(3);
phi = phi0 + XX(4);
kapa = kapa0 + XX(5);

% convert radian to degree
omega = omega*180/pi();
phi = phi*180/pi();
kapa = kapa*180/pi();


Orient = [omega; phi; kapa; By; Bz];