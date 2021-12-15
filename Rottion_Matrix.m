function M=Rottion_Matrix(omega,phi,kappa,units)

% Returns the rotation matrix given the three angles
% units: 0 rad (default)
%        1 grad
%        2 deg
%
%  Author: Ilias Kalisperakis 
%  2002
format long

if nargin==3 
    units=0;
end

if units ==1
    omega=omega*pi/200;
    phi=phi*pi/200;
    kappa=kappa*pi/200;
    
elseif units == 2
    omega=omega*pi()/180;
    phi=phi*pi()/180;
    kappa=kappa*pi()/180;
    
end

 Mx = [1, 0, 0;...
     0, cos(omega), sin(omega);...
     0, -sin(omega), cos(omega)];
 
 My = [cos(phi), 0, -sin(phi);...
     0, 1, 0;...
     sin(phi), 0, cos(phi)];

 Mz = [cos(kappa), sin(kappa), 0;...
     -sin(kappa), cos(kappa), 0;...
     0, 0, 1];
 
 M = Mz*My*Mx;
%  M1=[  (cos(kappa)*cos(phi))             (cos(kappa)*sin(phi)*sin(omega)+sin(kappa)*cos(omega))            (((-cos(kappa))*sin(phi)*cos(omega))+(sin(kappa)*sin(omega)))
%     (-sin(kappa)*cos(phi))      (((-sin(kappa))*sin(phi)*sin(omega))+(cos(kappa)*cos(omega)))                   (sin(kappa)*sin(phi)*cos(omega)+cos(kappa)*sin(omega))
%                   sin(phi)                                             (-cos(phi))*sin(omega)                                                      cos(phi)*cos(omega)];
 
end