% FUNCTION CORNERING: Finds the indices of the vertex of the
% cube of grid points surrounding the point (xx,yy,zz).
% ip = index of minimum x grid coord (min. x coord. = grid(1,ip)).
% jp = index of minimum y grid coord (min. y coord. = grid(2,jp)).
% kp = index of minimum z grid coord (min. z coord. = grid(3,kp)).

function [ip,jp,kp,flag]=cornering(xx,yy,zz,gridD)
    
    x(1)=xx;
    x(2)=yy;
    x(3)=zz;
    ngrid=[nnz(gridD(1,:)) nnz(gridD(2,:)) nnz(gridD(3,:))];
    flag=0;
    
    %If statement to check if the point is outside the grid
    if zz < gridD(3,1)
        [~,ip]=min(abs(xx-gridD(1,1:ngrid(1))));
        [~,jp]=min(abs(yy-gridD(2,1:ngrid(2))));
        kp=1;
        flag=1;
        no = zz;
        warning('Point is at depth=%f, too shallow',no);
        return
    elseif zz > gridD(3,ngrid(3))
        [~,ip]=min(abs(xx-gridD(1,1:ngrid(1))));
        [~,jp]=min(abs(yy-gridD(2,1:ngrid(2))));
        kp=ngrid(3);
        no = zz;
        warning('Point is at depth=%f, too deep',no);
        flag=2;
        return
    end
    if yy < gridD(2,1)
        [~,ip]=min(abs(xx-gridD(1,1:ngrid(1))));
        [~,kp]=min(abs(zz-gridD(3,1:ngrid(3))));
        jp=1;
        no = yy;
        warning('Point is at y=%f',no);
        flag=1;
        return
    elseif yy > gridD(2,ngrid(2))
        [~,ip]=min(abs(xx-gridD(1,1:ngrid(1))));
        [~,kp]=min(abs(zz-gridD(3,1:ngrid(3))));
        no = yy;
        warning('Point is at y=%f',no);
        jp=ngrid(2);
        flag=1;
        return
    end
    if xx < gridD(1,1)
        [~,jp]=min(abs(yy-gridD(2,1:ngrid(2))));
        [~,kp]=min(abs(zz-gridD(3,1:ngrid(3))));
        ip=1;
        no = xx;
        warning('Point is at x=%f',no);
        flag=1;
        return
    elseif xx > gridD(1,ngrid(1))
        [~,jp]=min(abs(yy-gridD(2,1:ngrid(2))));
        [~,kp]=min(abs(zz-gridD(3,1:ngrid(3))));
        ip=ngrid(1);
        no = xx;
        warning('Point is at x=%f',no);
        flag=1;
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
end
