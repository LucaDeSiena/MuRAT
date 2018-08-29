function v=vel(xx,yy,zz,ngrid,grid,pvel)
%     This routine finds the velocity at the point (xx,yy,zz) by
%     linear interpolation.
x(1)=xx;
x(2)=yy;
x(3)=zz;

%corner
if (xx < grid(1,1)) || (xx > grid(1,ngrid(1))) || ...
    (yy < grid(2,1)) || (yy > grid(2,ngrid(2))) || ...
    (zz < grid(3,1)) || (zz > grid(3,ngrid(3)))
    no = [xx,yy,zz];
    display(no);
    error('RAY BEND ERR: point is outside grid, called by corner');
end

% Loop over dimensions
indx=zeros(3,1);
for j=1:3
    indx(j)=ngrid(j)-1;
    while x(j)<grid(j,indx(j))
        indx(j)=indx(j)-1;
    end
end
ip=indx(1);
jp=indx(2);
kp=indx(3);

%vel
ip1=ip+1;
jp1=jp+1;
kp1=kp+1;
xd=grid(1,ip1)-grid(1,ip);
yd=grid(2,jp1)-grid(2,jp);     
zd=grid(3,kp1)-grid(3,kp);
xf=(xx-grid(1,ip))/xd;
yf=(yy-grid(2,jp))/yd;
zf=(zz-grid(3,kp))/zd;
v1=pvel(ip,jp,kp)+(pvel(ip1,jp,kp)-pvel(ip,jp,kp))*xf;
v2=pvel(ip,jp1,kp)+(pvel(ip1,jp1,kp)-pvel(ip,jp1,kp))*xf;
v3=v1+(v2-v1)*yf;
v4=pvel(ip,jp,kp1)+(pvel(ip1,jp,kp1)-pvel(ip,jp,kp1))*xf;
v5=pvel(ip,jp1,kp1)+(pvel(ip1,jp1,kp1)-pvel(ip,jp1,kp1))*xf;
v6=v4+(v5-v4)*yf;
v=v3+(v6-v3)*zf;     
