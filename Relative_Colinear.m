% Bismillah
function [Orient, XYZ]=Relative_Colinear(xy1, xy2, Base, xo , yo , f )
% about function: this function is used to do compute relative orientation parameters
% based on co-linear equation.
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
X0 = Base;
Y0 = 0;
Z0 = 0;
kapa0 = 0;
omega0 = 0;
phi0 = 0;

% compute parallax x
px = x1-x2;
% compute primery z.
n = size(x1,1);
m = n;
% ZA0 = Height - Base*f*ones(n,1)./px;
ZA0 =  - Base*f*ones(n,1)./px;
XA0 = ZA0.*(x1-xo)/(-f);
YA0 = ZA0.*(y1-yo)/(-f);
XYZ = [XA0, YA0, ZA0];

for i=1:40
    
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
     
     roundF2_omega = (f*ones(m,1)./(q.^2)).*(r.*(-M(3,3)*deltaY+M(3,2)*deltaZ) - q.*(-M(1,3)*deltaY+M(1,2)*deltaZ));
     roundF2_phi = (f*ones(m,1)./(q.^2)).*(r.*(Cph*deltaX + So*Sph*deltaY - Co*Sph*deltaZ)-...
                    q.*(-Sph*Ck*deltaX + So*Cph*Ck*deltaY -Co*Cph*Ck*deltaZ ));
     roundF2_kapa = (-f*ones(m,1)./q).*(M(2,1)*deltaX + M(2,2)*deltaY + M(2,3)*deltaZ);
     
     roundF2_Y0  = (-f*ones(m,1)./(q.^2)).*(r*M(3,2) - q*M(1,2));
     roundF2_Z0  = (-f*ones(m,1)./(q.^2)).*(r*M(3,3) - q*M(1,3));
     %
     roundG2_omega = (f*ones(m,1)./(q.^2)).*(s.*(-M(3,3)*deltaY+M(3,2)*deltaZ) - q.*(-M(2,3)*deltaY+M(2,2)*deltaZ));
     roundG2_phi = (f*ones(m,1)./(q.^2)).*(s.*(Cph*deltaX + So*Sph*deltaY - Co*Sph*deltaZ)-...
                    q.*(Sph*Sk*deltaX - So*Cph*Sk*deltaY + Co*Cph*Sk*deltaZ ));
     roundG2_kapa = (f*ones(m,1)./q).*(M(1,1)*deltaX + M(1,2)*deltaY + M(1,3)*deltaZ);
     
     roundG2_Y0  = (-f*ones(m,1)./(q.^2)).*(s*M(3,2) - q*M(2,2));
     roundG2_Z0  = (-f*ones(m,1)./(q.^2)).*(s*M(3,3) - q*M(2,3));
     % for 3d coordinates in model space
     roundF2_XA  = -(-f*ones(m,1)./(q.^2)).*(r*M(3,1) - q*M(1,1));
     roundF2_YA = -roundF2_Y0;
     roundF2_ZA = -roundF2_Z0;
     
     roundG2_XA  = -(-f*ones(m,1)./(q.^2)).*(s*M(3,1) - q*M(2,1));
     roundG2_YA = -roundG2_Y0;
     roundG2_ZA = -roundG2_Z0;

     F2X0 = xo -f * r./q;
     G2X0 = yo -f * s./q;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % for tie points im left image
     F1X0 = xo -f * XA0./ZA0;
     G1X0 = yo -f * YA0./ZA0;
     roundF1_XA  = -f*ones(m,1)./ZA0;
     roundF1_YA  = zeros(m,1);
     roundF1_ZA  = f*XA0./(ZA0.^2);
     
     roundG1_XA  = zeros(m,1);
     roundG1_YA  = -f*ones(m,1)./ZA0;
     roundG1_ZA  = f*YA0./(ZA0.^2);
     %
     A = zeros(4*n, 5+3*n);
     A(3:4:end,1:5)=[roundF2_omega, roundF2_phi, roundF2_kapa,  roundF2_Y0, roundF2_Z0];
     A(4:4:end,1:5)=[roundG2_omega, roundG2_phi, roundG2_kapa,  roundG2_Y0, roundG2_Z0];
     %
    for j=1:n
        % for tie points in left image
         A(4*j-3, 5+3*j-2) = roundF1_XA(j);
         A(4*j-3, 5+3*j-1) = roundF1_YA(j);
         A(4*j-3, 5+3*j) = roundF1_ZA(j);
         
         A(4*j-2, 5+3*j-2) = roundG1_XA(j);
         A(4*j-2, 5+3*j-1) = roundG1_YA(j);
         A(4*j-2, 5+3*j) = roundG1_ZA(j);
        % for tie points in right image
         A(4*j-1, 5+3*j-2) = roundF2_XA(j);
         A(4*j-1, 5+3*j-1) = roundF2_YA(j);
         A(4*j-1, 5+3*j) = roundF2_ZA(j);
         
         A(4*j, 5+3*j-2) = roundG2_XA(j);
         A(4*j, 5+3*j-1) = roundG2_YA(j);
         A(4*j, 5+3*j) = roundG2_ZA(j);
    end
     %
     L = zeros(4*n, 1);
     
     L(1:4:end) = x1 - F1X0 + roundF1_XA.*XA0 + roundF1_YA.*YA0 + roundF1_ZA.*ZA0;
     L(2:4:end) = y1 - G1X0 + roundG1_XA.*XA0 + roundG1_YA.*YA0 + roundG1_ZA.*ZA0;
     
     L(3:4:end) = x2 - F2X0 + roundF2_omega*omega0 + roundF2_phi*phi0 +...
                  roundF2_kapa*kapa0 + roundF2_Y0*Y0 + roundF2_Z0*Z0 +...
                  roundF2_XA.*XA0 + roundF2_YA.*YA0 + roundF2_ZA.*ZA0;
     L(4:4:end) = y2 - G2X0 + roundG2_omega*omega0 + roundG2_phi*phi0 +...
                  roundG2_kapa*kapa0 + roundG2_Y0*Y0 + roundG2_Z0*Z0 + ...
                  roundG2_XA.*XA0 + roundG2_YA.*YA0 + roundG2_ZA.*ZA0;
     %
     XX = inv(A'*A)*A'*L;%
   dx = XX(1:5,1)-[omega0; phi0; kapa0; Y0; Z0];
   dx1 = [XX(6:3:end), XX(7:3:end), XX(8:3:end)]-XYZ;
   dx2 = max(abs(dx1(:)));
   dx = max(abs(dx));
   
   if dx <0.0000001
       break;
   else
       omega0 = XX(1);
       phi0 = XX(2);
       kapa0 = XX(3);
       Y0 = XX(4);
       Z0 = XX(5);
       XA0 = XX(6:3:end,1);
       YA0 = XX(7:3:end,1);
       ZA0 = XX(8:3:end,1);
       XYZ = [XA0, YA0, ZA0];
   end
              
     
end
omega = XX(1)*180/pi();
phi = XX(2)*180/pi();
kapa = XX(3)*180/pi();
YC = XX(4);
ZC = XX(5);

XA0 = XX(6:3:end,1);
YA0 = XX(7:3:end,1);
ZA0 = XX(8:3:end,1);
XYZ = [XA0, YA0, ZA0];

Orient = [omega; phi; kapa; Base; YC; ZC];