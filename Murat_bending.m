function [xtemp,ytemp,ztemp,v]  =...
    Murat_bending(xtemp,ytemp,ztemp,gridD,v,pvel)
%BENDS the initial segment along the normal to the
% ray path tangent at each point by an optimal distance r.

xfac                            =   1;
n                               =   length(xtemp);

for k=2:(n-1)
    kk                          =   k-1;
    kkk                         =   k+1;
    
    % compute the normal direction of maximum gradient of velocity
    dx                          =   xtemp(kkk)-xtemp(kk);
    dy                          =   ytemp(kkk)-ytemp(kk);
    dz                          =   ztemp(kkk)-ztemp(kk);
    dn                          =   dx^2+dy^2+dz^2;
    ddn                         =   sqrt(dn);
    rdx                         =   dx/ddn;
    rdy                         =   dy/ddn;
    rdz                         =   dz/ddn;
    xk                          =   0.5*dx+xtemp(kk);
    yk                          =   0.5*dy+ytemp(kk);
    zk                          =   0.5*dz+ztemp(kk);
    
    vk                          =   Murat_velocity(xk,yk,zk,gridD,pvel);
    [vx,vy,vz]                  =   Murat_veldicity(xk,yk,zk,gridD,pvel);
    vrd                         =   vx*rdx+vy*rdy+vz*rdz;
    rvx                         =   vx-vrd*rdx;
    rvy                         =   vy-vrd*rdy;
    rvz                         =   vz-vrd*rdz;
    rvs                         =   sqrt(rvx*rvx+rvy*rvy+rvz*rvz);
    
    %  check for zero gradient
    if (rvs == 0)
        %  zero gradient - straight path
        xtemp(k)                =   xk;
        ytemp(k)                =   yk;
        ztemp(k)                =   zk;
        vk                      =   Murat_velocity(xk,yk,zk,gridD,pvel);
        v(k)                    =   vk;
    else
        rvx                     =   rvx/rvs;
        rvy                     =   rvy/rvs;
        rvz                     =   rvz/rvs;
        
        % compute the optimal distance r
        c                       =   (1/v(kkk)+1/v(kk))/2;
        c1                      =   (c*vk+1)/(4*c*rvs);
        rtemp                   =   -c1+sqrt(c1^2+dn/(8*c*vk));
        
        % compute the new points and distance of perturbations
        xxk                     =   xk+rvx*rtemp;
        yyk                     =   yk+rvy*rtemp;
        zzk                     =   zk+rvz*rtemp;
        
        %  convergence enhancement
        xxk                     =   xfac*(xxk-xtemp(k))+xtemp(k);
        yyk                     =   xfac*(yyk-ytemp(k))+ytemp(k);
        zzk                     =   xfac*(zzk-ztemp(k))+ztemp(k);
        
        xtemp(k)                =   xxk;
        ytemp(k)                =   yyk;
        ztemp(k)                =   zzk;
        v(k)                    =   Murat_velocity(xxk,yyk,zzk,gridD,pvel);
    end
end
