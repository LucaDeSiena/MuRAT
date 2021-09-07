function [vx,vy,vz]             =   Murat_veldicity(xx,yy,zz,gridD,pvel)
% function [vx,vy,vz]             =   Murat_veldicity(xx,yy,zz,gridD,pvel)
%
% CALCULATES the gradient of the velocity in x, y, and z directions with linear interpolation
%
% Input parameters:
%    xx:        x point
%    yy:        y point
%    zz:        z point
%    gridD:     grid of ray tracing
%    pvel:      velocity model for ray tracing
%
% Output parameters:
%    vx:        dv/dx
%    vy:        dv/dx
%    vz:        dv/dx

xGrid                       =   gridD.x;
yGrid                       =   gridD.y;
zGrid                       =   gridD.z;

[ip,jp,kp,flag]             =   Murat_cornering(xx,yy,zz,gridD);

if flag>0
    vx                      =   0;
    vy                      =   0;
    vz                      =   0;
    return
end

ip1                         =   ip+1;
jp1                         =   jp+1;
kp1                         =   kp+1;
xd                          =   xGrid(ip1) - xGrid(ip);
yd                          =   yGrid(jp1) - yGrid(jp);
zd                          =   zGrid(kp1) - zGrid(kp);

xf                          =   (xx - xGrid(ip))/xd;
yf                          =   (yy - yGrid(jp))/yd;
zf                          =   (zz - zGrid(kp))/zd;

xf1                         =   1 - xf;
yf1                         =   1 - yf;
zf1                         =   1 - zf;

% calculate derivatives of velocity
vx                          =(yf1*zf1*(pvel(ip1,jp,kp)-pvel(ip,jp,kp))+...
    yf*zf1*(pvel(ip1,jp1,kp)-pvel(ip,jp1,kp))+ ...
    yf1*zf*(pvel(ip1,jp,kp1)-pvel(ip,jp,kp1))+ ...
    yf*zf*(pvel(ip1,jp1,kp1)-pvel(ip,jp1,kp1)))/xd;

vy                          =(xf1*zf1*(pvel(ip,jp1,kp)-pvel(ip,jp,kp))+...
    xf*zf1*(pvel(ip1,jp1,kp)-pvel(ip1,jp,kp))+ ...
    xf1*zf*(pvel(ip,jp1,kp1)-pvel(ip,jp,kp1))+ ...
    xf*zf*(pvel(ip1,jp1,kp1)-pvel(ip1,jp,kp1)))/yd;

vz                          =(xf1*yf1*(pvel(ip,jp,kp1)-pvel(ip,jp,kp))+...
    xf*yf1*(pvel(ip1,jp,kp1)-pvel(ip1,jp,kp))+ ...
    xf1*yf*(pvel(ip,jp1,kp1)-pvel(ip,jp1,kp))+ ...
    xf*yf*(pvel(ip1,jp1,kp1)-pvel(ip1,jp1,kp)))/zd;
end