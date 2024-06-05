function [modvProp,modvP,modvEqSpace,modvPlot]...
                            =...
    Murat_modv3D(FPath,FLabel,modvOriginal,origin,mLat,mLon,nLat,nLon,nzc)
% function [modvP,modvPropagation,modvEqS,modvPlot]	=...
%     Murat_modv3D(FPath,FLabel,modvOriginal,origin,mLat,mLon,nLat,nLon)
% CREATES the velocity models when using a 3D as input
%
% Input parameters:
%    FPath,FLabel:          Paths to save test velocity figure
%    modvOriginal:          original 3D velocity model
%    origin:                origin of the spatial grid
%    mLat:                  min/max latitudes of earthquakes
%    mLon:                  min/max longitudes of earthquakes
%    nLat:                  number nodes latitude
%    nLon:                  number nodes longitude
%    nzc:                   number nodes vertical
%
% Output parameters:
%    modvP:                 propagation velocity model in standard format
%    modvPropagation:       propagation velocity model in matrix format
%    modvEqS:               equally-spaced inversion velocity model
%    modvPlot:              equally-spaced ploting velocity model (lat/lon)

%% Rescale original lat/lon modV in cartesian coords
% Transform nodes velocity model in meters with WGS84
wgs84                       =   wgs84Ellipsoid("m");
[d_o,az_o]                  =   distance(origin(1),origin(2),...
    modvOriginal(:,1),modvOriginal(:,2),wgs84);

% Transform min/max lat/lon earthquakes in meters with WGS84
[d_ll,az_ll]                =   distance(origin(1),origin(2),...
    mLat',mLon',wgs84);
mLatLon                     =   [d_ll.*sin(az_ll*2*pi/360)...
    d_ll.*cos(az_ll*2*pi/360)];

%% Original velocity model in meters - uneven spacing
modvOrC                     =   sortrows([d_o.*sin(az_o*2*pi/360)...
    d_o.*cos(az_o*2*pi/360) modvOriginal(:,3:4)],[1,2,-3]);
xOrC                        =   modvOrC(:,1);
yOrC                        =   modvOrC(:,2);
zOrC                        =   modvOrC(:,3);
vOrC                        =   modvOrC(:,4);


%% Velocity model for ray tracing - inversion
% For ray-tracing, you need to space evenly with the same number of
% x-y nodes as in input - here we extend to min and max quake locations
lModvP                      = [min(mLatLon(1,1),min(xOrC))...
    max(mLatLon(2,1),max(xOrC)); min(mLatLon(1,2),min(yOrC))...
    max(mLatLon(2,2),max(yOrC))];

xM1                         =   linspace(lModvP(1,1),lModvP(1,2),nLon)';
yM1                         =   linspace(lModvP(2,1),lModvP(2,2),nLat)';
zM1                         =...
    sort(linspace(min(zOrC),max(zOrC),nzc),'descend');


[XEqS,YEqS,ZEqS,modvEqSpace,modvEqS] =...
                            Murat_rescale(xOrC,yOrC,zOrC,vOrC,xM1,yM1,zM1);

% Change to lat/lon for plotting
d                           =...
    sqrt(modvEqSpace(:,1).^2 + modvEqSpace(:,2).^2);
az                          =...
    atan(modvEqSpace(:,1)./modvEqSpace(:,2))*360/2/pi;
[lat2,lon2]                 =   reckon(origin(1),origin(2),d,az,wgs84);
modvPlot                    =   [lon2 lat2 modvEqSpace(:,3:4)];

%Interpolated model for ray-tracing - half grid step of original
resol2x                     =   abs(xM1(2)-xM1(1))/2;
resol2y                     =   abs(yM1(2)-yM1(1))/2;
resol2z                     =   abs(zM1(2)-zM1(1))/2;

xp                          =   xM1(1):resol2x:xM1(end);
yp                          =   yM1(1):resol2y:yM1(end);
zp                          =   zM1(1):-resol2z:zM1(end);

[Xq,Yq,Zq,modvProp,modvP]   =...
                            Murat_rescale(xOrC,yOrC,zOrC,vOrC,xp,yp,zp);

Murat_test_vel_models(FPath,FLabel,modvOriginal,modvOrC,...
        XEqS,YEqS,ZEqS,modvEqS,Xq,Yq,Zq,modvP);
end