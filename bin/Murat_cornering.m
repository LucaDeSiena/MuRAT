function [ip,jp,kp,flag]    =   Murat_cornering(xx,yy,zz,gridD)
% function [ip,jp,kp,flag]    =   Murat_cornering(xx,yy,zz,gridD)
%
% CHECKS that the point is in the grid and
% FINDS the indices of the vertex of the cube [ip,jp,kp] of grid points
%   surrounding the point [xx,yy,zz].
%
% Input parameters:
%    xx:        x coord of point West;
%    yy:        y coord of point North;
%    zz:        z coord of point ;
%    gridD:     propagation grid;
%
% Output parameters:
%    ip:        index of minimum x grid coord - West;
%    jp:        index of minimum y grid coord - South;
%    kp:        index of maximum z grid coord - shallowest.
%    flag:      flags if ray outside of grid

% Coordinates of the propagation grid are unwrapped
xGrid   =   gridD.x;
yGrid   =   gridD.y;
zGrid   =   gridD.z;

% If point exceeds grid extrema
flag    =   0;

% If statement to check if the point is outside the grid
if zz < min(zGrid)
    [~,ip]  =   min(abs(xx-xGrid));
    [~,jp]  =   min(abs(yy-yGrid));
    [~,kp]  =   min(zGrid);
    flag    =   1;
    return
elseif zz > max(zGrid)
    [~,ip]  =   min(abs(xx-xGrid));
    [~,jp]  =   min(abs(yy-yGrid));
    [~,kp]  =   max(zGrid);
    flag    =   1;
    return
end
if yy < min(yGrid)
    [~,ip]  =   min(abs(xx-xGrid));
    [~,jp]  =   min(yGrid);
    [~,kp]  =   min(abs(zz-zGrid));
    flag    =   1;
    return
elseif yy > max(yGrid)
    [~,ip]  =   min(abs(xx-xGrid));
    [~,jp]  =   max(yGrid);
    [~,kp]  =   min(abs(zz-zGrid));
    flag    =   1;
    return
end
if xx < min(xGrid)
    [~,ip]  =   min(xGrid);
    [~,jp]  =   min(abs(yy-yGrid));
    [~,kp]  =   min(abs(zz-zGrid));
    flag=1;
    return
elseif xx > max(xGrid)
    [~,ip]  =   max(xGrid);
    [~,jp]  =   min(abs(yy-yGrid));
    [~,kp]  =   min(abs(zz-zGrid));
    flag    =   1;
    return
end

% Finds index - remember convention on shallowest SW point!
% A point on a grid node (including the west, south and top edges) belongs
%   to the cell that starts there; the last node maps to the last cell.
ip          =   min(find(xGrid<=xx,1,'last'),numel(xGrid)-1);
jp          =   min(find(yGrid<=yy,1,'last'),numel(yGrid)-1);
kp          =   min(find(zGrid>=zz,1,'last'),numel(zGrid)-1);

end
