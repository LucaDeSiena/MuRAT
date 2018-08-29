function [ta,tra]=travel(xtemp,ytemp,ztemp,ngrid,grid,v,pvel)

n = length(xtemp);
xd=xtemp(2)-xtemp(1);
yd=ytemp(2)-ytemp(1);
zd=ztemp(2)-ztemp(1);
ds=sqrt(xd^2+yd^2+zd^2);
if ds <= 5000
    tra=ds*((1/v(1)+1/v(2))/2);
    ta(2)=tra;
    for i=3:n
        i1=i-1;
        xd=xtemp(i)-xtemp(i1);
        yd=ytemp(i)-ytemp(i1);
        zd=ztemp(i)-ztemp(i1);
        ds=sqrt(xd^2+yd^2+zd^2);
        tra=tra+ds/((v(i)+v(i1))/2);
        ta(i)=tra;
    end
else
    tra=0;
    for i=2:n
        i1=i-1;
        xd=xtemp(i)-xtemp(i1);
        yd=ytemp(i)-ytemp(i1);
        zd=ztemp(i)-ztemp(i1);
        ds=sqrt(xd^2+yd^2+zd^2);
        xcos=xd/ds;
        ycos=yd/ds;
        zcos=zd/ds;
        is=fix(ds/5000);
        ds=ds/is;
        for j=1:is
            xx=xtemp(i1)+(j-0.5)*ds*xcos;
            yy=ytemp(i1)+(j-0.5)*ds*ycos;
            zz=ztemp(i1)+(j-0.5)*ds*zcos;
            vp=vel(xx,yy,zz,ngrid,grid,pvel);
            tra=tra+ds/vp;
        end
        ta(i)=tra;
    end
end