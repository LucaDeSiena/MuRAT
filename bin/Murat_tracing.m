function ma     =   Murat_tracing(ray,gridP,pvel)
% function ma   =   Murat_tracing(ray,gridP,pvel)
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
maxit       =   100; 

% Max bend points, 1 direction 
maxpoints   =   10000; 

% Three dimensions
max3        =   3*maxpoints; 

% Nodes of the grid
xNodes      =   gridP.x;
yNodes      =   gridP.y;
zNodes      =   gridP.z;

% Max distance one cell
dconv       =   sqrt((xNodes(2)-xNodes(1))^2 +...
    (yNodes(2)-yNodes(1))^2 + (zNodes(2)-zNodes(1))^2);

% Time limit, in seconds, pvel is in km
tlim        =   dconv/mean(mean(mean(pvel*1000))); 

% Start bending procedure from the extrema of the ray segment
xtemp       =   [ray(1,1) (ray(1,1)+ray(1,2))/2 ray(1,2)];
ytemp       =   [ray(2,1) (ray(2,1)+ray(2,2))/2 ray(2,2)];
ztemp       =   [ray(3,1) (ray(3,1)+ray(3,2))/2 ray(3,2)];

n           =   3;
v           =   zeros(n,1);

for i=1:n
    v(i)    =   Murat_velocity(xtemp(i), ytemp(i), ztemp(i), gridP, pvel);
end

% Computes travel time
[ta,tra]    =   Murat_traveling(xtemp,ytemp,ztemp,gridP,v,pvel);
ttsave      =   0;

for it = 1:maxit
    prevT   =   tra;
    
    [xtemp,ytemp,ztemp,v] =  Murat_bending(xtemp,ytemp,ztemp,gridP,v,pvel);
    
    [ta,tra]    =   Murat_traveling(xtemp,ytemp,ztemp,gridP,v,pvel);
    
    deltat      =   prevT-tra;
    
    if deltat < tlim
        % Check convergence and distance criterium
        dtemp   =   sqrt((xtemp(2)-xtemp(1))^2 + (ytemp(2)-ytemp(1))^2 ...
            + (ztemp(2)-ztemp(1))^2);
        
        if abs(ttsave-tra) < tlim && dtemp < dconv/10
            break
        
        else
            % As the travel time decrease is smaller than travel-time
            %   treshold (tlim) and/or the travel time increased, the
            %   code doubles the number of segments:
            ttsave  =   tra;
            nnew    =   2*n-1;
            if(nnew > max3)
                error('Too many points, reduce them or increase max3');
            end
            
            % allocate new arrays
            xt      =   zeros(1,nnew);
            yt      =   zeros(1,nnew);
            zt      =   zeros(1,nnew);

            % place original points into odd indices (1,3,5,...)
            oddIdx  =   1:2:nnew;
            xt(oddIdx)  =   xtemp;
            yt(oddIdx)  =   ytemp;
            zt(oddIdx)  =   ztemp;

            % fill midpoints into even indices
            evenIdx     =   2:2:(nnew-1);
            xt(evenIdx) =   0.5 * (xt(oddIdx(1:end-1)) + xt(oddIdx(2:end)));
            yt(evenIdx) =   0.5 * (yt(oddIdx(1:end-1)) + yt(oddIdx(2:end)));
            zt(evenIdx) =   0.5 * (zt(oddIdx(1:end-1)) + zt(oddIdx(2:end)));

            xtemp   = xt; ytemp = yt; ztemp = zt;
            n       = nnew;

            % recompute velocities at all nodes
            v = zeros(n,1);
            for i = 1:n
                v(i)=Murat_velocity(xtemp(i),ytemp(i),ztemp(i),gridP,pvel);
            end
            % reset ta so next Murat_traveling recalculates correctly
            ta      =   zeros(1,n);
        end
    end
end
% assemble output: columns [index, x, y, z, travel_time_seconds]
ma      = zeros(n,5);
ma(:,1) = (1:n)';
ma(:,2) = xtemp(:);
ma(:,3) = ytemp(:);
ma(:,4) = ztemp(:);
ma(:,5) = ta(:)' / 1000;   % preserve original division by 1000 (units)
end