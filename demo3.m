% bismillah
clc
clear
A= zeros(11,5);
A(1,1)=1;
A(2,3)=1;
A(3,4)=1;
A(4,1)=-1; A(4,3)=1;
A(5,1)=-1; A(5,2)=1;
A(6,1)=-1; A(6,4)=1;
A(7,2)=-1; A(7,3)=1;
A(8,2)=-1;
A(9,2)=-1; A(9,5)=1;
A(10,5)=-1;
A(11,4)=1; A(11,5)=-1;
%
L=[102.644; 99.849; 101.516; -2.803; 4.691; -1.135; -7.478;...
    -107.337; -3.799; -103.534; -2.033];
x= inv(A'*A)*A'*L;
v=A*x-L;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
xA=1000;
yA=1000;

xB=1200;
yB=990;

xC=1280;
yC=1180;

x0=1100;
y0=1150;

L1 = 180.287;
L2 = 188.668;
L3 = 182.482;

for i=1:20
    f1x0 = sqrt((xA-x0)^2+(yA-y0)^2);
    f2x0 = sqrt((xB-x0)^2+(yB-y0)^2);
    f3x0 = sqrt((xC-x0)^2+(yC-y0)^2);
    %
    roun_f1x1 = -(xA-x0)/L1;
    roun_f1x2 = -(yA-y0)/L1;
    
    roun_f2x1 = -(xB-x0)/L2;
    roun_f2x2 = -(yB-y0)/L2;
    
    roun_f3x1 = -(xC-x0)/L3;
    roun_f3x2 = -(yC-y0)/L3;
    %
    A = [roun_f1x1 roun_f1x2;
         roun_f2x1 roun_f2x2;
         roun_f3x1 roun_f3x2];
    %
    L=[L1-f1x0 + roun_f1x1*x0 + roun_f1x2*y0;
       L2-f2x0 + roun_f2x1*x0 + roun_f2x2*y0;
       L3-f3x0 + roun_f3x1*x0 + roun_f3x2*y0];
   x = inv(A'*A)*A'*L;
   dx = max(abs(x(1)-x0), abs(x(2)-y0));
   
   if dx <0.000001
       break;
   else
       x0 = x(1);
       y0 = x(2);
   end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
% assignment 3: linear modelling
A= [2,1
    3,1
    6,1
    11,1
    13,1
    17,1
    19,1
    25,1];
y=[-1.932
    0.076
    6.074
    16.039
    20.066
    28.017
    32.071
    44.003];
x= inv(A'*A)*A'*y;
    


