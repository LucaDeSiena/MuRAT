function [vx,vy,vz]     =   Murat_veldicity(xx,yy,zz,gridD,pvel)
% function [vx,vy,vz]   =   Murat_veldicity(xx,yy,zz,gridD,pvel)
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

xg          =   gridD.x;
yg          =   gridD.y;
zg          =   gridD.z;

% Outside the grid the velocity is constant along the normal to the nearest
%   grid face (see Murat_velocity): use the gradient on that face and set
%   the components along the directions where the point is outside to 0.
outX        =   xx < min(xg) || xx > max(xg);
outY        =   yy < min(yg) || yy > max(yg);
outZ        =   zz < min(zg) || zz > max(zg);
xx          =   min(max(xx,min(xg)),max(xg));
yy          =   min(max(yy,min(yg)),max(yg));
zz          =   min(max(zz,min(zg)),max(zg));

[ip,jp,kp]  =   Murat_cornering(xx,yy,zz,gridD);

ip1         =   ip+1;
jp1         =   jp+1;
kp1         =   kp+1;

xd          =   xg(ip1) - xg(ip);
yd          =   yg(jp1) - yg(jp);
zd          =   zg(kp1) - zg(kp);

xf          =   (xx - xg(ip))/xd;
yf          =   (yy - yg(jp))/yd;
zf          =   (zz - zg(kp))/zd;

xf1         =   1 - xf;
yf1         =   1 - yf;
zf1         =   1 - zf;

% corner values (consistent with original pvel(jp,ip,kp) ordering)
v000        =   pvel(jp,  ip,  kp);
v100        =   pvel(jp,  ip1, kp);
v010        =   pvel(jp1, ip,  kp);
v110        =   pvel(jp1, ip1, kp);
v001        =   pvel(jp,  ip,  kp1);
v101        =   pvel(jp,  ip1, kp1);
v011        =   pvel(jp1, ip,  kp1);
v111        =   pvel(jp1, ip1, kp1);

% dv/dx: differences along x between corresponding corners, weighted by y/z
vx = (yf1*zf1*(v100 - v000) + ... % lower z, lower y
    yf*zf1 *(v110 - v010) + ...   % lower z, upper y
    yf1*zf  *(v101 - v001) + ...  % upper z, lower y
    yf*zf  *(v111 - v011)) / xd;  % upper z, upper y

% dv/dy: differences along y between corresponding corners, weighted by x/z
vy = (xf1*zf1*(v010 - v000) + ... % lower z, lower x
    xf*zf1 *(v110 - v100) + ...   % lower z, upper x
    xf1*zf  *(v011 - v001) + ...  % upper z, lower x
    xf*zf  *(v111 - v101)) / yd;  % upper z, upper x

% dv/dz: differences along z between corresponding corners, weighted by x/y
vz = (xf1*yf1*(v001 - v000) + ... % lower y, lower x
    xf*yf1 *(v101 - v100) + ...   % lower y, upper x
    xf1*yf  *(v011 - v010) + ...  % upper y, lower x
    xf*yf  *(v111 - v110)) / zd;  % upper y, upper x

if outX; vx = 0; end
if outY; vy = 0; end
if outZ; vz = 0; end

end