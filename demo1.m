% Bismillah
clc
clear

f = 0.152;% meters
tilt = 5;% degrees
r = 0.110; % meters
landa = 30;% degrees
dr = r^2 * sind(tilt) * cosd(landa)/(f - r * sind(tilt) * cosd(landa))
dr_mm = dr *1000
