%FINDS the indices of the vertex of the cube of grid points surrounding the
% point [xx,yy,zz].
% ip = index of minimum x grid coord - West;
% jp = index of minimum y grid coord - South;
% kp = index of maximum z grid coord - shallowest.

function [ip,jp,kp,flag]    =   Murat_cornering(xx,yy,zz,gridD)

% Coordinates of the propagation grid are unwrapped
xGrid                       =   gridD.x;
yGrid                       =   gridD.y;
zGrid                       =   gridD.z;

% If point exceeds grid extrema
flag                        =   0;

%If statement to check if the point is outside the gridg
if zz < min(zGrid)
    [~,ip]                  =   min(abs(xx-xGrid));
    [~,jp]                  =   min(abs(yy-yGrid));
    [~,kp]                  =   min(zGrid);
    flag                    =   1;
    %     no = zz;
    %     warning('Point is at depth=%f, too deep',no);
    return
elseif zz > max(zGrid)
    [~,ip]                  =   min(abs(xx-xGrid));
    [~,jp]                  =   min(abs(yy-yGrid));
    [~,kp]                  =   max(zGrid);
    flag                    =   1;
    %     no = zz;
    %     warning('Point is at depth=%f, too shallow',no);
    return
end
if yy < min(yGrid)
    [~,ip]                  =   min(abs(xx-xGrid));
    [~,jp]                  =   min(yGrid);
    [~,kp]                  =   min(abs(zz-zGrid));
    flag                    =   1;
    %     no = yy;
    %     warning('Point is at y=%f',no);
    return
elseif yy > max(yGrid)
    [~,ip]                  =   min(abs(xx-xGrid));
    [~,jp]                  =   max(yGrid);
    [~,kp]                  =   min(abs(zz-zGrid));
    flag                    =   1;
    %     no = yy;
    %     warning('Point is at y=%f',no);
    return
end
if xx < min(xGrid)
    [~,ip]                  =   min(xGrid);
    [~,jp]                  =   min(abs(yy-yGrid));
    [~,kp]                  =   min(abs(zz-zGrid));
    flag=1;
    %     no = xx;
    %     warning('Point is at x=%f',no);
    return
elseif xx > max(xGrid)
    [~,ip]                  =   max(xGrid);
    [~,jp]                  =   min(abs(yy-yGrid));
    [~,kp]                  =   min(abs(zz-zGrid));
    flag                    =   1;
    %     no = xx;
    %     warning('Point is at x=%f',no);
    return
end

% Finds index - remember convention on shallowest SW point!
ip                          =   find(xGrid<xx,1,'last');
jp                          =   find(yGrid<yy,1,'last');
kp                          =   find(zGrid>zz,1,'last');

end
