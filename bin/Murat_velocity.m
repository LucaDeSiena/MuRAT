% FUNCTION Murat_velocity: It finds the velocity at the point (xx,yy,zz) by
% linear interpolation.

function v  =   Murat_velocity(xx,yy,zz,gridD,pvel)
%
% CALCULATES the velocity at xx, yy, and zz by linear interpolation
%
% Input parameters:
%    xx:        x point
%    yy:        y point
%    zz:        z point
%    gridD:     grid of ray tracing
%    pvel:      velocity model for ray tracing
%
% Output parameters:
%    v:         velocity

xGrid       =   gridD.x;
yGrid       =   gridD.y;
zGrid       =   gridD.z;

[ip,jp,kp,flag]	=   Murat_cornering(xx,yy,zz,gridD);
if flag>0
    v       =   pvel(jp,ip,kp);
    return
end

ip1         =   ip+1;
jp1         =   jp+1;
kp1         =   kp+1;

xd          =   xGrid(ip1) - xGrid(ip);
yd          =   yGrid(jp1) - yGrid(jp);
zd          =   zGrid(kp1) - zGrid(kp);

xf          =   (xx - xGrid(ip))/xd;
yf          =   (yy - yGrid(jp))/yd;
zf          =   (zz - zGrid(kp))/zd;

% corner values (consistent indexing: pvel(j,y,x) as original)
v000        =   pvel(jp,  ip,  kp);
v100        =   pvel(jp,  ip1, kp);
v010        =   pvel(jp1, ip,  kp);
v110        =   pvel(jp1, ip1, kp);
v001        =   pvel(jp,  ip,  kp1);
v101        =   pvel(jp,  ip1, kp1);
v011        =   pvel(jp1, ip,  kp1);
v111        =   pvel(jp1, ip1, kp1);

% trilinear interpolation
c00         =   v000*(1 - xf) + v100*xf;
c10         =   v010*(1 - xf) + v110*xf;
c01         =   v001*(1 - xf) + v101*xf;
c11         =   v011*(1 - xf) + v111*xf;

c0          =   c00*(1 - yf) + c10*yf;
c1          =    c01*(1 - yf) + c11*yf;

v           =   c0*(1 - zf) + c1*zf;
end