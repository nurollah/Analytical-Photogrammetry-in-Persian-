% Bismillah
clc
clear

% question 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part A:
% inputs:
f = 0.152;% % focal length( meters)
tilt = 5; % degrees
tilt = tilt * pi()/180;% radian
r = 0.11; % raduis in meters
landa = 30;% degrees
landa = landa * pi()/180;% radian

% process
dr_tilt = (r^2 * sin(tilt) * cos(landa))/(f - r * sin(tilt) * cos(landa)); % tilt displacement (meters)
dr_tilt = dr_tilt*1000;% tilt displacement (milimeters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part B:
h = 15;% height (meters)
r = 0.1; % raduis ( meters)
H = 800; % flying height (meters)
dr_Height = r*h/H; % height displacement (meters)
dr_Height = dr_Height *1000; % height displacement (milimeters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part C:
H = 2000;% flying height ( meters)
t = 1/500;% shutter speed (second)
ha= 0 ; % height ( meters)
f = 0.180;% focal length( meters)
V = 800 *1000/3600;% Velocity (m/s)
dr_IM = V * t * f /(H-ha); % image motion displacement (meters)
dr_IM = dr_IM *10^6; % image motion displacement (micron)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part D:
f = 0.180;% focal length( meters)
r = 0.06; % raduis ( meters)
H = 1500;% flying height ( meters)
R = 6372230;% raduis for earth ( meters)
dr_refraction = H * r^3 / (2 * R * f^2); % Refraction displacement (meters)
dr_refraction = dr_refraction * 10^6; % Refraction displacement (micron)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part E:
x_pp = 0.008;% the x coordinate of principal point (mm)
y_pp = -0.001;% the y coordinate of principal point (mm)
x_old = 62.579;% the x coordinate of interest point (mm)
y_old = -80.916;% the y coordinate of interest point (mm)
k1 = 0.2296;% the first coefficient
k2 = -35.89;% the second coefficient
k3 = 1018;% the third coefficient
k4 = 12.1;% the fourth coefficient
rmm = sqrt((x_old - x_pp)^2 + (y_old - y_pp)^2 );% raduis ( milimeters)
r = rmm/1000;% raduis ( meters)
dr = k1*r + k2*r^3 + k3*r^5 + k4*r^7;% radial displacement (mm)

x_new = x_old * (1-dr/rmm);% the correcrted coordinate (mm)
y_new = y_old * (1-dr/rmm);% the correcrted coordinate (mm)

