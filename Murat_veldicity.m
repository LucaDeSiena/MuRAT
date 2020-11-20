% FUNCTION VELDICITY: gradient of the velocity in x, y, and z directions at
% (xx,yy,zz): vx=dv/dx; vy=dv/dy; vz=dv/dz. It uses linear interpolation
% for speed.
function [vx,vy,vz]=Murat_veldicity(xx,yy,zz,gridD,pvel)
    
    [ip,jp,kp,flag]= Murat_cornering(xx,yy,zz,gridD);
    if flag>0
        vx=0;
        vy=0;
        vz=0;
        return
    end
    
    ip1=ip+1;
    jp1=jp+1;
    kp1=kp+1;
    xd=gridD(1,ip1)-gridD(1,ip);
    yd=gridD(2,jp1)-gridD(2,jp);
    zd=gridD(3,kp1)-gridD(3,kp);
    xf=(xx-gridD(1,ip))/xd;
    yf=(yy-gridD(2,jp))/yd;
    zf=(zz-gridD(3,kp))/zd;
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