function DD_coord               =   Murat_DDcoordinates(origin,ending,nLat,nLong,nzc)
%%
% function DD_coord               =   Murat_DDcoordinates(origin,ending,nLat,nLong,nzc)
%
% Recalculate grid points coordinates in decimal degrees starting from 
% origin and adding the mean spacing between points (in decimal degrees).
%
% Input parameters:
%    origin:   vector of coordinates of grid origin point [lon,lat,depth]
%    ending:   vector of coordinates of grid ending point [lon,lat,depth]
%    nLat:     number of grid nodes in latitude
%    nLong:    number of grid nodes in longitude
%    nzc:      number of grid nodes in depth
%
% Output parameters:
%    DD_coord: matrix of coordinates [x,y,z] for each node in the Murat
%    order

%% 
% coordinate vectors
yLat        = linspace(origin(1), ending(1), nLat);   % latitude (j)
xLong       = linspace(origin(2), ending(2), nLong);  % longitude (i)
zDep        = linspace(origin(3), ending(3), nzc);    % depth (k)

% produce grids with same indexing as original nested loops:
[X, Y, Z]   = ndgrid(xLong, yLat, zDep);

% Generating matrix of coordinates
DD_coord    = [X(:), Y(:), Z(:)];

end