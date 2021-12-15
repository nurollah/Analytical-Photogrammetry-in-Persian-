% Bismillah
function [Par, residuals]=Affine_3D_2D(XYZ, xy)
format long;
m = size(XYZ,1);
A = zeros(2*m, 8);
L = zeros(2*m, 1);
%
L(1:2:end,1) = xy(:,1);
L(2:2:end,1) = xy(:,2);

A(1:2:end,1) = XYZ(:,1);
A(1:2:end,2) = XYZ(:,2);
A(1:2:end,3) = XYZ(:,3);
A(1:2:end,4) = 1;
A(2:2:end,5) = XYZ(:,1);
A(2:2:end,6) = XYZ(:,2);
A(2:2:end,7) = XYZ(:,3);
A(2:2:end,8) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Par = pinv(A'*A)*A'*L;
residuals1 = A*Par - L;
residuals = [residuals1(1:2:end), residuals1(2:2:end)];
