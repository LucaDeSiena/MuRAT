function [xtemp,ytemp,ztemp,v]  =...
    Murat_bending(xtemp,ytemp,ztemp,gridD,v,pvel)
% function [xtemp,ytemp,ztemp,v]=...
%  Murat_bending(xtemp,ytemp,ztemp,gridD,v,pvel)
% 
% BENDS the initial segment along the normal to the
%   ray path tangent at each point by an optimal distance r.
%   Uses the standard ray-bending approach - Block, 1991.
%
% Input parameters:
%    xtemp:     x before bending
%    ytemp:     y before bending
%    ztemp:     z before bending
%    gridD:     grid of ray tracing
%    v:         initial velocity of the ray
%    pvel:      velocity model for ray tracing
%
% Output parameters:
%    xtemp:     x after bending
%    ytemp:     y after bending
%    ztemp:     z after bending
%    v:         final velocity of the ray

xfac                =   1;
n                   =   length(xtemp);

for k=2:(n-1)
    kk              =   k-1;
    kkk             =   k+1;
    
    % compute the normal direction of maximum gradient of velocity
    % segment vector and its length
    dx              =   xtemp(kkk) - xtemp(kk);
    dy              =   ytemp(kkk) - ytemp(kk);
    dz              =   ztemp(kkk) - ztemp(kk);
    dn              =   dx*dx + dy*dy + dz*dz;
    ddn             =   sqrt(dn);
    inv_ddn         =   1 / ddn;
    rdx             =   dx * inv_ddn;
    rdy             =   dy * inv_ddn;
    rdz             =   dz * inv_ddn;

    % midpoint
    xk              =   xtemp(kk) + 0.5*dx;
    yk              =   ytemp(kk) + 0.5*dy;
    zk              =   ztemp(kk) + 0.5*dz;
    
    % velocity and gradient (vel dicity) at midpoint
    vk              =   Murat_velocity(xk,yk,zk,gridD,pvel);
    [vx,vy,vz]      =   Murat_veldicity(xk,yk,zk,gridD,pvel);
    
    % gradient component perpendicular to tangent
    vrd             =   vx*rdx + vy*rdy + vz*rdz;
    rvx             =   vx - vrd*rdx;
    rvy             =   vy - vrd*rdy;
    rvz             =   vz - vrd*rdz;
    rvs             =   sqrt(rvx*rvx + rvy*rvy + rvz*rvz);
    
    %  check for zero gradient
    if (rvs == 0)
        %  zero gradient - straight path
        xtemp(k)    =   xk;
        ytemp(k)    =   yk;
        ztemp(k)    =   zk;
        v(k)        =   vk;
    else
        inv_rvs     =   1 / rvs;
        rvx         =   rvx * inv_rvs;
        rvy         =   rvy * inv_rvs;
        rvz         =   rvz * inv_rvs;
        
        % compute the optimal distance r
        c           =   (1/v(kkk) + 1/v(kk))*0.5;
        c1          =   (c*vk+1) / (4*c*rvs);
        rtemp       =   -c1 + sqrt(c1^2 + dn/(8*c*vk));
        
        % compute the new points and distance of perturbations
        xxk         =   xk + rvx*rtemp;
        yyk         =   yk + rvy*rtemp;
        zzk         =   zk + rvz*rtemp;
        
        %  convergence enhancement
        xxk         =   xfac * (xxk-xtemp(k)) + xtemp(k);
        yyk         =   xfac * (yyk-ytemp(k)) + ytemp(k);
        zzk         =   xfac * (zzk-ztemp(k)) + ztemp(k);
        
        xtemp(k)    =   xxk;
        ytemp(k)    =   yyk;
        ztemp(k)    =   zzk;
        v(k)        =   Murat_velocity(xxk,yyk,zzk,gridD,pvel);
    end
end