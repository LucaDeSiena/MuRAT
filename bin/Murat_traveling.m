function [ta,tra]       =   Murat_traveling(xtemp,ytemp,ztemp,gridD,v,pvel)
% function [ta,tra]   =   Murat_traveling(xtemp,ytemp,ztemp,gridD,v,pvel)
%
% COMPUTES travel time between bent points.
%
% Input parameters:
%    xtemp:     x before bending
%    ytemp:     y before bending
%    ztemp:     z before bending
%    gridD:     grid of ray tracing
%    v:         initial velocity of the ray
%    pvel:      velocity model for ray tracing in km
%
% Output parameters:
%    ta:        travel time index
%    tra:       travel time

% Coordinates of the propagation grid are unwrapped
xGrid                   =   gridD.x;
yGrid                   =   gridD.y;
zGrid                   =   gridD.z;

maxGrid                 =   max([xGrid(2)-xGrid(1)...
    yGrid(2)-yGrid(1) zGrid(2)-zGrid(1)]);

% Max travel distance limit
distlim                 =   5*(maxGrid);
pvel                    =   pvel*1000;

n                       =   length(xtemp);
xd                      =   xtemp(2)-xtemp(1);
yd                      =   ytemp(2)-ytemp(1);
zd                      =   ztemp(2)-ztemp(1);
ds                      =   sqrt(xd^2+yd^2+zd^2);
if ds <= distlim
    tra                 =   ds*((1/v(1)+1/v(2))/2);
    ta(2)               =   tra;
    for i = 3:n
        i1              =   i-1;
        xd              =   xtemp(i)-xtemp(i1);
        yd              =   ytemp(i)-ytemp(i1);
        zd              =   ztemp(i)-ztemp(i1);
        ds              =   sqrt(xd^2+yd^2+zd^2);
        tra             =   tra+ds/((v(i)+v(i1))/2);
        ta(i)           =   tra;
    end
else
    tra                 =   0;
    ta                  =   zeros(1,n);
    for i = 2:n
        i1              =   i-1;
        xd              =   xtemp(i)-xtemp(i1);
        yd              =   ytemp(i)-ytemp(i1);
        zd              =   ztemp(i)-ztemp(i1);
        ds              =   sqrt(xd^2+yd^2+zd^2);
        xcos            =   xd/ds;
        ycos            =   yd/ds;
        zcos            =   zd/ds;
        is              =   fix(ds/distlim);
        ds              =   ds/is;
        for j = 1:is
            xx          =   xtemp(i1)+(j-0.5)*ds*xcos;
            yy          =   ytemp(i1)+(j-0.5)*ds*ycos;
            zz          =   ztemp(i1)+(j-0.5)*ds*zcos;
            vp          =   Murat_velocity(xx,yy,zz,gridD,pvel);
            tra         =   tra+ds/vp;
        end
        ta(i)           =   tra;
    end
end
end