function [vx, vy, vz, t] = extractLOH1data(station)

Nt = 450; 
% ns = 9; number of stations
dt = 0.02;

t = (0:Nt-1)*dt;

L = load('LOH1data.mat');

vx = L.x((station-1)*Nt+1:station*Nt);
vy = L.y((station-1)*Nt+1:station*Nt);
vz = L.z((station-1)*Nt+1:station*Nt);


vx = deconv_conv(vx,t);
vy = deconv_conv(vy,t);  
vz = deconv_conv(vz,t);
end