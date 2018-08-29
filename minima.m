function ma=minima(ray,ngrid,grid,pvel)

%  This routine finds the minimum path using pseudo-bending
%  Initial path is straight ray, defined by two endpoints.
%  Straight ray mid-point is used as the third point neccessary
%  for ray bending. 
maxit=1000;
tlim = 0.1;
blkpar=2000;
blkpar3=3*blkpar;
dconv=20;
%  next line for more exact tracing
tlim2=tlim;
% BENT RAY

n=3;

xtemp(1)=ray(1,1);
ytemp(1)=ray(2,1);
ztemp(1)=ray(3,1);
xtemp(3)=ray(1,2);
ytemp(3)=ray(2,2);
ztemp(3)=ray(3,2);
xtemp(2)=(xtemp(1)+xtemp(3))/2;
ytemp(2)=(ytemp(1)+ytemp(3))/2;
ztemp(2)=(ztemp(1)+ztemp(3))/2;

for i=1:n
    xp=xtemp(i);
    yp=ytemp(i);
    zp=ztemp(i);
    vp=vel(xp,yp,zp,ngrid,grid,pvel);
    v(i,1)=vp;
end
[ta,tra]=travel(xtemp,ytemp,ztemp,ngrid,grid,v,pvel);
ttsave =0;
for j=1:maxit
    taa=tra;
    [xtemp,ytemp,ztemp,v]=bend(xtemp,ytemp,ztemp,ngrid,grid,v,pvel);
    [ta,tra]=travel(xtemp,ytemp,ztemp,ngrid,grid,v,pvel);
    deltat=taa-tra;
    if deltat < tlim2
        % check for convergence
        dtemp=sqrt((xtemp(2)-xtemp(1))^2+(ytemp(2)-ytemp(1))^2 ...
            +(ztemp(2)-ztemp(1))^2);
        if abs(ttsave-tra) < tlim && dtemp < dconv
            break
        else
            ttsave=tra;
            % If the travel time decrease was smaller than the travel-time
            % improvement parameter (tlim) (or if the travel time increased),
            % double the number of path segments:
            ntemp=2*n-1;
            if(ntemp > blkpar3)
               error('ERROR: minima6 - number of points gt blkpar3');
            end
            ta(1)=0;
            for i=1:n
                l=n-i+1;
                k=2*l-1;
                xtemp(k)=xtemp(l);
                ytemp(k)=ytemp(l);
                ztemp(k)=ztemp(l);
                ta(k)=ta(l);
            end
            for i=2:2:(2*n-2)
                xtemp(i)=(xtemp(i-1)+xtemp(i+1))/2;
                ytemp(i)=(ytemp(i-1)+ytemp(i+1))/2;
                ztemp(i)=(ztemp(i-1)+ztemp(i+1))/2;
            end
            n=2*n-1;
            for i=1:n
                xp=xtemp(i);
                yp=ytemp(i);
                zp=ztemp(i);
                vp=vel(xp,yp,zp,ngrid,grid,pvel);
                v(i)=vp;
            end
        end
    end
end
for i=1:n
    x(i)=xtemp(i);
    y(i)=ytemp(i);
    z(i)=ztemp(i);
    ma(i,:) = [i,x(i),y(i),z(i),ta(i)];
end
end