% FUNCTION VELOCITY: It finds the velocity at the point (xx,yy,zz) by
% linear interpolation.

function v=velocity(xx,yy,zz,gridD,pvel)
    x(1)=xx;
    x(2)=yy;
    x(3)=zz;
    ngrid=[nnz(gridD(1,:)) nnz(gridD(2,:)) nnz(gridD(3,:))];
    
    % Finding the corner of the block
    [ip,jp,kp,flag]=cornering(xx,yy,zz,gridD);
    if flag>0
        v=pvel(ip,jp,kp);
        return
    end
    
    % Loop over dimensions
    indx=zeros(3,1);
    for j=1:3
        indx(j)=ngrid(j)-1;
        while x(j)<gridD(j,indx(j))
            indx(j)=indx(j)-1;
        end
    end
    ip=indx(1);
    jp=indx(2);
    kp=indx(3);
    
    ip1=ip+1;
    jp1=jp+1;
    kp1=kp+1;
    xd=gridD(1,ip1)-gridD(1,ip);
    yd=gridD(2,jp1)-gridD(2,jp);
    zd=gridD(3,kp1)-gridD(3,kp);
    xf=(xx-gridD(1,ip))/xd;
    yf=(yy-gridD(2,jp))/yd;
    zf=(zz-gridD(3,kp))/zd;
    v1=pvel(ip,jp,kp)+(pvel(ip1,jp,kp)-pvel(ip,jp,kp))*xf;
    v2=pvel(ip,jp1,kp)+(pvel(ip1,jp1,kp)-pvel(ip,jp1,kp))*xf;
    v3=v1+(v2-v1)*yf;
    v4=pvel(ip,jp,kp1)+(pvel(ip1,jp,kp1)-pvel(ip,jp,kp1))*xf;
    v5=pvel(ip,jp1,kp1)+(pvel(ip1,jp1,kp1)-pvel(ip,jp1,kp1))*xf;
    v6=v4+(v5-v4)*yf;
    v=v3+(v6-v3)*zf;
end
