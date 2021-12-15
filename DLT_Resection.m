% Bismillah
function [Par, residuals]=DLT_Resection(XYZ, xy)
format long;
m = size(XYZ,1);
A = zeros(2*m, 11);
L = zeros(2*m, 1);
%
L(1:2:end,1) = xy(:,1);
L(2:2:end,1) = xy(:,2);

A(1:2:end,1) = XYZ(:,1);
A(1:2:end,2) = XYZ(:,2);
A(1:2:end,3) = XYZ(:,3);
A(1:2:end,4) = 1;
%
A(2:2:end,5) = XYZ(:,1);
A(2:2:end,6) = XYZ(:,2);
A(2:2:end,7) = XYZ(:,3);
A(2:2:end,8) = 1;
%
A(1:2:end,9) = -xy(:,1) .* XYZ(:,1);% -x*X
A(1:2:end,10) = -xy(:,1) .* XYZ(:,2);% -x*Y
A(1:2:end,11) = -xy(:,1) .* XYZ(:,3);% -x*Z
%
A(2:2:end,9) = -xy(:,2) .* XYZ(:,1);% -y*X
A(2:2:end,10) = -xy(:,2) .* XYZ(:,2);% -y*Y
A(2:2:end,11) = -xy(:,2) .* XYZ(:,3);% -y*Z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Par = pinv(A'*A)*A'*L;
residuals1 = A*Par - L;
residuals = [residuals1(1:2:end), residuals1(2:2:end)];
