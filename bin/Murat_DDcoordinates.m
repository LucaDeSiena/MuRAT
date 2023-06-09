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
% Calculating increments
lat_incr                        =   (ending(1)-origin(1))/(nLat-1);
lon_incr                        =   (ending(2)-origin(2))/(nLong-1);
dep_incr                        =   (ending(3)-origin(3))/(nzc-1);
% Creating vectors for x,y,z
xLong                           =   zeros(nLong,1);
yLat                            =   zeros(nLat,1);
zDep                            =   zeros(nzc,1);
lon=origin(2);
for i=1:nLong
    xLong(i)=lon;
    lon=lon+lon_incr;
end
lat=origin(1);
for i=1:nLat
    yLat(i)=lat;
    lat=lat+lat_incr;
end
dep=origin(3);
for i=1:nzc
    zDep(i)=dep;
    dep=dep+dep_incr;
end
% Generating matrix of coordinates
DD_coord=zeros(nLong*nLat*nzc,3);
index=0;
for i=1:nLong
    for j=1:nLat
        for k=1:nzc
            index       =   index+1;                
            DD_coord(index,1:3)   =   [xLong(i) yLat(j) zDep(k)];
        end
    end
end
end