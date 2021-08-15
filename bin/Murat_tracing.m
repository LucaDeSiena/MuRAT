function ma                 =   Murat_tracing(ray,gridP,pvel)
% function ma                 =   Murat_tracing(ray,gridP,pvel)
% 
% TRACES the minimum path between source and receiver using pseudo-bending
%
% Input parameters:
%    ray:       input ray
%    gridP:     grid for ray tracing
%    gridD:     grid of ray tracing
%    pvel:      velocity model for ray tracing
%
% Output parameters:
%    ma:        ray in Murat format
% 
% Structure:
% The initial path is a straight segment. The segment mid-point is used for
%   bending. The procedure repeats for a maximum number of points, or until
%   the changes are below a certain treshold. 

% Maximum number of iterations allowed
maxit                       =   100; 

% Max bend points, 1 direction 
maxpoints                   =   10000; 

% Three dimensions
max3                        =   3*maxpoints; 

% Nodes of the grid
xNodes                      =   gridP.x;
yNodes                      =   gridP.y;
zNodes                      =   gridP.z;

% Max distance one cell
dconv                       =   sqrt((xNodes(2)-xNodes(1))^2 +...
    (yNodes(2)-yNodes(1))^2 + (zNodes(2)-zNodes(1))^2);

% Time limit, in seconds, pvel is in km
tlim                        =   dconv/mean(mean(mean(pvel*1000))); 

% Start bending procedure from the extrema of the ray segment
xtemp                       =   [ray(1,1) (ray(1,1)+ray(1,2))/2 ray(1,2)];
ytemp                       =   [ray(2,1) (ray(2,1)+ray(2,2))/2 ray(2,2)];
ztemp                       =   [ray(3,1) (ray(3,1)+ray(3,2))/2 ray(3,2)];

n                           =   3;
v                           =   zeros(n,1);

for i=1:n
    xp                      =   xtemp(i);
    yp                      =   ytemp(i);
    zp                      =   ztemp(i);
    vp                      =   Murat_velocity(xp,yp,zp,gridP,pvel);
    v(i,1)                  =   vp;
end

% Computes travel time
[ta,tra]                    =...
    Murat_traveling(xtemp,ytemp,ztemp,gridP,v,pvel);

ttsave                      =   0;
for j=1:maxit
    taa                     =   tra;
    
    [xtemp,ytemp,ztemp,v]   =...
        Murat_bending(xtemp,ytemp,ztemp,gridP,v,pvel);
    
    [ta,tra]                =...
        Murat_traveling(xtemp,ytemp,ztemp,gridP,v,pvel);
    
    deltat                  =   taa-tra;
    
    % Check for convergence
    if deltat < tlim
        dtemp               =...
            sqrt((xtemp(2)-xtemp(1))^2 + (ytemp(2)-ytemp(1))^2 ...
            + (ztemp(2)-ztemp(1))^2);
        
        if abs(ttsave-tra) < tlim && dtemp < dconv
            break
        
        else
            % As the travel time decrease is smaller than travel-time
            %   treshold (tlim) and/or the travel time increased, the
            %   code doubles the number of segments:
            ttsave          =   tra;
            ntemp           =   2*n-1;
            if(ntemp > max3)
                error('Too many points, reduce them or increase max3');
            
            end
            
            ta(1)=0;
            for i=1:n
                l           =   n-i+1;
                k           =   2*l-1;
                xtemp(k)    =   xtemp(l);
                ytemp(k)    =   ytemp(l);
                ztemp(k)    =   ztemp(l);
                ta(k)       =   ta(l);
            end
            
            for i=2:2:(2*n-2)
                xtemp(i)    =   (xtemp(i-1)+xtemp(i+1))/2;
                ytemp(i)    =   (ytemp(i-1)+ytemp(i+1))/2;
                ztemp(i)    =   (ztemp(i-1)+ztemp(i+1))/2;
            end
            
            n               =   2*n-1;
            for i=1:n
                xp          =   xtemp(i);
                yp          =   ytemp(i);
                zp          =   ztemp(i);
                vp          =   Murat_velocity(xp,yp,zp,gridP,pvel);
                v(i)        =   vp;
            end
        end
    end
end

% Final ray, including travel time info
ma(1:n,1:5)                 =   [(1:n)',xtemp',ytemp',ztemp',ta'/1000]; 
end