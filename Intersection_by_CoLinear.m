% Bismillah
function [XYZ]=Intersection_by_CoLinear(xy1, xy2, exterior, interior)
%
xo = interior(1,1);
yo = interior(1,2);
f = interior(1,3);
%
omega1=exterior(1,1); 
phi1 = exterior(1,2); 
kappa1 = exterior(1,3);
X01 = exterior(1,4);
Y01 = exterior(1,5);
Z01 = exterior(1,6);
%
omega2 =exterior(2,1); 
phi2 = exterior(2,2); 
kappa2 = exterior(2,3);
X02 = exterior(2,4);
Y02 = exterior(2,5);
Z02 = exterior(2,6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M1 = Rottion_Matrix(omega1, phi1, kappa1, 2);
M2 = Rottion_Matrix(omega2, phi2, kappa2, 2);
%
R = M1';
%
x1 = xy1(1,1);
y1 = xy1(1,2);
x2 = xy2(1,1);
y2 = xy2(1,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute initial values
Base = sqrt((X02 - X01)^2+(Y02 - Y01)^2);
px = x1 - x2;% parallax

ZA0 = (Z01 + Z02)/2 - Base*f/px;% initial value for Z

XA0 = X01 + (ZA0 -Z01)*(R(1,1)*(x1-xo) + R(1,2)*(y1-yo) + R(1,3)*(-f))/...
                        (R(3,1)*(x1-xo) + R(3,2)*(y1-yo) + R(3,3)*(-f));
                    
YA0 = Y01 + (ZA0 -Z01)*(R(2,1)*(x1-xo) + R(2,2)*(y1-yo) + R(2,3)*(-f))/...
                        (R(3,1)*(x1-xo) + R(3,2)*(y1-yo) + R(3,3)*(-f)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:20
    F1X0 = xo -f * (M1(1,1)*(XA0-X01) + M1(1,2)*(YA0-Y01) + M1(1,3)*(ZA0-Z01))/...
                    (M1(3,1)*(XA0-X01) + M1(3,2)*(YA0-Y01) + M1(3,3)*(ZA0-Z01));
    G1X0 = yo -f * (M1(2,1)*(XA0-X01) + M1(2,2)*(YA0-Y01) + M1(2,3)*(ZA0-Z01))/...
                    (M1(3,1)*(XA0-X01) + M1(3,2)*(YA0-Y01) + M1(3,3)*(ZA0-Z01));
    F2X0 = xo -f * (M2(1,1)*(XA0-X02) + M2(1,2)*(YA0-Y02) + M2(1,3)*(ZA0-Z02))/...
                    (M2(3,1)*(XA0-X02) + M2(3,2)*(YA0-Y02) + M2(3,3)*(ZA0-Z02));
    G2X0 = yo -f * (M2(2,1)*(XA0-X02) + M2(2,2)*(YA0-Y02) + M2(2,3)*(ZA0-Z02))/...
                    (M2(3,1)*(XA0-X02) + M2(3,2)*(YA0-Y02) + M2(3,3)*(ZA0-Z02));
    %
    round_F1X = (-M1(3,1)*(x1-xo)-f*M1(1,1))/...
                (M1(3,1)*(XA0-X01) + M1(3,2)*(YA0-Y01) + M1(3,3)*(ZA0-Z01));
    round_F1Y = (-M1(3,2)*(x1-xo)-f*M1(1,2))/...
                (M1(3,1)*(XA0-X01) + M1(3,2)*(YA0-Y01) + M1(3,3)*(ZA0-Z01));
    round_F1Z = (-M1(3,3)*(x1-xo)-f*M1(1,3))/...
                (M1(3,1)*(XA0-X01) + M1(3,2)*(YA0-Y01) + M1(3,3)*(ZA0-Z01));
    
    round_G1X = (-M1(3,1)*(y1-yo)-f*M1(2,1))/...
                (M1(3,1)*(XA0-X01) + M1(3,2)*(YA0-Y01) + M1(3,3)*(ZA0-Z01));
    round_G1Y = (-M1(3,2)*(y1-yo)-f*M1(2,2))/...
                (M1(3,1)*(XA0-X01) + M1(3,2)*(YA0-Y01) + M1(3,3)*(ZA0-Z01));
    round_G1Z = (-M1(3,3)*(y1-yo)-f*M1(2,3))/...
                (M1(3,1)*(XA0-X01) + M1(3,2)*(YA0-Y01) + M1(3,3)*(ZA0-Z01));
    
    round_F2X = (-M2(3,1)*(x2-xo)-f*M2(1,1))/...
                (M2(3,1)*(XA0-X02) + M2(3,2)*(YA0-Y02) + M2(3,3)*(ZA0-Z02));
    round_F2Y = (-M2(3,2)*(x2-xo)-f*M2(1,2))/...
                (M2(3,1)*(XA0-X02) + M2(3,2)*(YA0-Y02) + M2(3,3)*(ZA0-Z02));
    round_F2Z = (-M2(3,3)*(x2-xo)-f*M2(1,3))/...
                (M2(3,1)*(XA0-X02) + M2(3,2)*(YA0-Y02) + M2(3,3)*(ZA0-Z02));
    
    round_G2X = (-M2(3,1)*(y2-yo)-f*M2(2,1))/...
                (M2(3,1)*(XA0-X02) + M2(3,2)*(YA0-Y02) + M2(3,3)*(ZA0-Z02));
    round_G2Y = (-M2(3,2)*(y2-yo)-f*M2(2,2))/...
                (M2(3,1)*(XA0-X02) + M2(3,2)*(YA0-Y02) + M2(3,3)*(ZA0-Z02));
    round_G2Z = (-M2(3,3)*(y2-yo)-f*M2(2,3))/...
                (M2(3,1)*(XA0-X02) + M2(3,2)*(YA0-Y02) + M2(3,3)*(ZA0-Z02));
    
    %
    A = [round_F1X round_F1Y round_F1Z;
         round_G1X round_G1Y round_G1Z;
         round_F2X round_F2Y round_F2Z;
         round_G2X round_G2Y round_G2Z];
    %
    L=[x1-F1X0 + round_F1X*XA0 + round_F1Y*YA0 + round_F1Z*ZA0;
       y1-G1X0 + round_G1X*XA0 + round_G1Y*YA0 + round_G1Z*ZA0;
       x2-F2X0 + round_F2X*XA0 + round_F2Y*YA0 + round_F2Z*ZA0;
       y2-G2X0 + round_G2X*XA0 + round_G2Y*YA0 + round_G2Z*ZA0];
   XX = inv(A'*A)*A'*L;
   dx = max(abs(XX(1)-XA0), max( abs(XX(2)-YA0), abs(XX(3)-ZA0)));
   
   if dx <0.000001
       break;
   else
       XA0 = XX(1);
       YA0 = XX(2);
       ZA0 = XX(3);
   end
    
end

XYZ =XX';

