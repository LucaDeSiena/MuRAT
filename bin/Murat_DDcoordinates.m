%% Recalculate grid point in decimal degrees
function DD_coord               =   Murat_DDcoordinates(Murat)
%%
% Importing all the necessary inputs
origin                          =   Murat.input.origin;
ending                          =   Murat.input.end;
nLat                            =   Murat.input.gridLat;
nLong                           =   Murat.input.gridLong;
nDep                            =   Murat.input.gridZ;
% Calculating increments
lat_incr                        =   (ending(1)-origin(1))/(nLat-1);
lon_incr                        =   (ending(2)-origin(2))/(nLong-1);
dep_incr                        =   (ending(3)-origin(3))/(nDep-1);
% Creating vectors for x,y,z
xLong                           =   zeros(nLong,1);
yLat                            =   zeros(nLat,1);
zDep                            =   zeros(nDep,1);
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
for i=1:nDep
    zDep(i)=dep;
    dep=dep+dep_incr;
end
% Generating matrix of coordinates
DD_coord=zeros(nLong*nLat*nDep,3);
index=0;
for i=1:nLong
    for j=1:nLat
        for k=1:nDep
            index       =   index+1;                
            DD_coord(index,1:3)   =   [xLong(i) yLat(j) zDep(k)];
        end
    end
end
end