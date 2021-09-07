% FUNCTION Murat_velocity: It finds the velocity at the point (xx,yy,zz) by
% linear interpolation.

function v      =   Murat_velocity(xx,yy,zz,gridD,pvel)

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

xGrid           =   gridD.x;
yGrid           =   gridD.y;
zGrid           =   gridD.z;

[ip,jp,kp,flag]	=   Murat_cornering(xx,yy,zz,gridD);
if flag>0
    v           =   pvel(ip,jp,kp);
    return
end

ip1             =   ip+1;
jp1             =   jp+1;
kp1             =   kp+1;

xd              =   xGrid(ip1) - xGrid(ip);
yd              =   yGrid(jp1) - yGrid(jp);
zd              =   zGrid(kp1) - zGrid(kp);

xf              =   (xx - xGrid(ip))/xd;
yf              =   (yy - yGrid(jp))/yd;
zf              =   (zz - zGrid(kp))/zd;

v1              =   pvel(ip,jp,kp) + (pvel(ip1,jp,kp) -pvel(ip,jp,kp))*xf;

v2              =   pvel(ip,jp1,kp)+(pvel(ip1,jp1,kp) -pvel(ip,jp1,kp))*xf;

v3              =   v1 + (v2 - v1)*yf;

v4              =   pvel(ip,jp,kp1)+(pvel(ip1,jp,kp1) -pvel(ip,jp,kp1))*xf;

v5              =   pvel(ip,jp1,kp1)+...
    (pvel(ip1,jp1,kp1) - pvel(ip,jp1,kp1))*xf;

v6              =   v4 + (v5 - v4)*yf;

v               =   v3 + (v6 - v3)*zf;
end