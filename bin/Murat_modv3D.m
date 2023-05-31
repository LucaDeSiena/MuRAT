function [modvP,modvI,modvIP,pvel]	=...
    Murat_modv3D(modvXYZ,modvOriginal,origin,flagTest)
% function [modvP,modvI,modvIP,pvel]	=...
%     Murat_modv3D(modvXYZ,modvOriginal,origin,flagTest)
% CREATES the velocity models when using a 3D as input
%
% Input parameters:
%    modvXYZ:               output grids of MuRAT_fold
%    modvOriginal:          original 3D velocity model
%    origin:                origin of the spatial grid
%    flagTest:              choose if testing the 3D with images
%
% Output parameters:
%    modvP:                 propagation velocity model
%    modvI:                 inversion velocity model
%    modvIP:                inversion velocity model in matrix
%    pvel:                  propagation velocity model in matrix

% Nodes of the inversion velocity model
xM                                  =   unique(modvXYZ(:,1));
yM                                  =   unique(modvXYZ(:,2));
zM                                  =   sort(unique(modvXYZ(:,3)),'descend');
[X,Y,Z]                             =   meshgrid(xM,yM,zM);

% modvDEG                             =...
%     [modvOriginal(:,1)-origin(1) modvOriginal(:,2)-origin(2)...
%     modvOriginal(:,3:4)];
% 
% % Switch as [lat,lon] is [y,x] - the velocity in Murat format
% modv_o                              =...
%     sortrows([deg2km(modvDEG(:,2))*1000 deg2km(modvDEG(:,1))*1000 ...
%     modvDEG(:,3:4)],[1,2,-3]);

% Transform nodes velocity model in meters with WGS84
wgs84                               =   wgs84Ellipsoid("m");
[d,az]                              =   distance(origin(1),origin(2),...
                                modvOriginal(:,1),modvOriginal(:,2),wgs84);

modv_o                              =   [d.*sin(az*2*pi/360)...
                                  d.*cos(az*2*pi/360) modvOriginal(:,3:4)];
[XM,YM,ZM,VM]                        =   Murat_fold(xM,yM,zM,modv_o(:,4));

% Space evenly with the same number of x-y nodes as in input
minModv                             =   min(modv_o);
maxModv                             =   max(modv_o);

xM1                                 =...
    linspace(minModv(2),maxModv(2),nLong)';
yM1                                 =...
    linspace(minModv(1),maxModv(1),nLat)';
zM1                                  =...
    sort(unique(modv_o(:,3)),'descend');

[XM1,YM1,ZM1]                        =   meshgrid(xM1,yM1,zM1);
modvIPM                              =   interp3(XM,YM,ZM,VM,XM1,YM1,ZM1);

if find(isnan(modvIPM))
    modvIPM                          =   inpaintn(modvIPM);
end
% EDIT
% Nodes of the original velocity model
xD                                  =   unique(modv_o(:,1));
yD                                  =   unique(modv_o(:,2));
zD                                  =...
    sort(unique(modv_o(:,3)),'descend');

% Uses meshgrid - fold switches again for interpolation and plot
[XD,YD,ZD,V]                        =   Murat_fold(xD,yD,zD,modv_o(:,4));
% Interpolated 3D inversion model
modvIP                              =   interp3(XD,YD,ZD,V,X,Y,Z);

%In case limits outside of the grid: interpolate
if find(isnan(modvIP))
    modvIP                          =   inpaintn(modvIP);
end


%Interpolated model for ray-tracing - half grid step of original
resol2x                             =   abs(xM(2)-xM(1))/2;
resol2y                             =   abs(yM(2)-yM(1))/2;
resol2z                             =   abs(zM(2)-zM(1))/2;

% Interpolated vectors
xp                                  =   xM(1):resol2x:xM(end);
yp                                  =   yM(1):resol2y:yM(end);
zp                                  =   zM(1):-resol2z:zM(end);

% Meshes for interpolation
[Xq,Yq,Zq]                          =   meshgrid(xp',yp',zp');

% pvel is the propagation velocity model in matrix format
pvel                                =   interp3(XD,YD,ZD,V,Xq,Yq,Zq);
if find(isnan(pvel))
    pvel                            =   inpaintn(pvel);
end

% Standard unpacking of the velocity model for propagation
modvP                               =   Murat_unfold(Xq,Yq,Zq,pvel);
modvI                               =   Murat_unfold(X,Y,Z,modvIP);
    
if isequal(flagTest,1)
    Murat_test_vel_models(modvOriginal,...
        modv_o,modvI,X,Y,Z,pvel,Xq,Yq,Zq,origin)
end

end