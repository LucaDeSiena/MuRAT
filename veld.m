function [vx,vy,vz]=veld(xx,yy,zz,ddn,ngrid,grid,pvel)
% This routine computes the derivatives of velocity
% in x, y, and z directions - uses linear interpolation
% (xx,yy,,zz) - point at which derivatives are to be computed
%  vx=dv/dx; vy=dv/dy; vz=dv/dz
delta=ddn/100;
if delta <= 25
    [ip,jp,kp]= corner(xx,yy,zz,ngrid,grid);

else
    vp=vel(xx-delta,yy-delta,zz-delta,ngrid,grid,pvel);
    veloc(1,1,1)=vp;
    vp=vel(xx+delta,yy-delta,zz-delta,ngrid,grid,pvel);
    veloc(2,1,1)=vp;
    vp = vel(xx-delta,yy+delta,zz-delta,ngrid,grid,pvel);
    veloc(1,2,1)=vp;
    vp = vel(xx+delta,yy+delta,zz-delta,ngrid,grid,pvel);
    veloc(2,2,1)=vp;
    vp = vel(xx-delta,yy-delta,zz+delta,ngrid,grid,pvel);
    veloc(1,1,2)=vp;
    vp = vel(xx+delta,yy-delta,zz+delta,ngrid,grid,pvel);
    veloc(2,1,2)=vp;
    vp = vel(xx-delta,yy+delta,zz+delta,ngrid,grid,pvel);
    veloc(1,2,2)=vp;
    vp = vel(xx+delta,yy+delta,zz+delta,ngrid,grid,pvel);
    veloc(2,2,2)=vp;
    % calculate derivatives of velocity
    vx=((veloc(2,1,1)-veloc(1,1,1))+(veloc(2,2,1)-veloc(1,2,1)) ...
        +(veloc(2,1,2)-veloc(1,1,2))+(veloc(2,2,2)-veloc(1,2,2)))/(8.*delta);

    vy=((veloc(1,2,1)-veloc(1,1,1))+(veloc(2,2,1)-veloc(2,1,1)) ...
        +(veloc(1,2,2)-veloc(1,1,2))+(veloc(2,2,2)-veloc(2,1,2)))/(8.*delta);
   
    vz=((veloc(1,1,2)-veloc(1,1,1))+(veloc(2,1,2)-veloc(2,1,1)) ...
        +(veloc(1,2,2)-veloc(1,2,1))+(veloc(2,2,2)-veloc(2,2,1)))/(8.*delta);

    [ip,jp,kp]= corner(xx,yy,zz,ngrid,grid);
end

ip1=ip+1;
jp1=jp+1;
kp1=kp+1;
xd=grid(1,ip1)-grid(1,ip);
yd=grid(2,jp1)-grid(2,jp);
zd=grid(3,kp1)-grid(3,kp);
xf=(xx-grid(1,ip))/xd;
yf=(yy-grid(2,jp))/yd;
zf=(zz-grid(3,kp))/zd;
xf1=1-xf;
yf1=1-yf;
zf1=1-zf;

% calculate derivatives of velocity
vx=(yf1*zf1*(pvel(ip1,jp,kp)-pvel(ip,jp,kp))+ ...
    yf*zf1*(pvel(ip1,jp1,kp)-pvel(ip,jp1,kp))+ ...
    yf1*zf*(pvel(ip1,jp,kp1)-pvel(ip,jp,kp1))+ ...
    yf*zf*(pvel(ip1,jp1,kp1)-pvel(ip,jp1,kp1)))/xd;

vy=(xf1*zf1*(pvel(ip,jp1,kp)-pvel(ip,jp,kp))+ ...
    xf*zf1*(pvel(ip1,jp1,kp)-pvel(ip1,jp,kp))+ ...
    xf1*zf*(pvel(ip,jp1,kp1)-pvel(ip,jp,kp1))+ ...
    xf*zf*(pvel(ip1,jp1,kp1)-pvel(ip1,jp,kp1)))/yd;

vz=(xf1*yf1*(pvel(ip,jp,kp1)-pvel(ip,jp,kp))+ ...
    xf*yf1*(pvel(ip1,jp,kp1)-pvel(ip1,jp,kp))+ ...
    xf1*yf*(pvel(ip,jp1,kp1)-pvel(ip,jp1,kp))+ ...
    xf*yf*(pvel(ip1,jp1,kp1)-pvel(ip1,jp1,kp)))/zd;
end