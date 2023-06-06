function testModv           =...
Murat_test_vel_models(FPath,FLabel,modvOriginal,modvOrCartesian,...
        XEqSpace,YEqSpace,ZEqSpace,modvEqS,Xq,Yq,Zq,modvP)
% function [modvP,modvI,modvIP,pvel]	=...
%     Murat_modv3D(modvXYZ,modvOriginal,origin,flagTest)
% TESTS correct recovery of velocity model when using a 3D as input
% Input parameters
%    modvOriginal:                  original velocity model
%    modvOrCartesian:               cartesian velocity model
%    XEqSpace,YEqSpace,ZEqSpace:    coords equally-spaced velocity model
%    modvEqS:                       equally-spaced velocity model
%    Xq,Yq,Zq:                      coords propagation velocity model
%    modvP:                         propagation velocity model
%
% Output figure:
%    testModv:          single figure showing original, cartesian,
%                       equally-spaced and propagation velocity model
    
testModv                    =   figure('Name','3DVelocityTest',...
    'NumberTitle','off','Position',[20,400,1200,1000],'visible','off');

% first look at original model
% use min and maximum lat / lon of model
lat_res                     =   unique(modvOriginal(:,1));
lon_res                     =   unique(modvOriginal(:,2));
depth_res                   =   unique(modvOriginal(:,3));

sliceLat                    =   mean(lat_res);
sliceLon                    =   mean(lon_res);
sliceDepth                  =   mean(depth_res);

%3D
[XOrModel,YOrModel,ZOrModel]=   meshgrid(lon_res,lat_res,depth_res);
VOrModel                    =   griddata(modvOriginal(:,2),...
    modvOriginal(:,1),modvOriginal(:,3),modvOriginal(:,4),...
    XOrModel,YOrModel,ZOrModel);

[color]                     =   inferno(100);

subplot(2,2,1)
slice(XOrModel,YOrModel,ZOrModel, VOrModel, sliceLon, sliceLat, sliceDepth)
colormap(color)
view(0,90)
hcb                         =   colorbar;
title(hcb,'V (km/s)')
title('Original velocity model')


%% check out centered model
% use min and maximum lon / lat of model
lon_res2                    =   unique(modvOrCartesian(:,1));
lat_res2                    =   unique(modvOrCartesian(:,2));
depth_res2                  =   unique(modvOrCartesian(:,3));

%3D
[Xmodv_o,Ymodv_o,Zmodv_o]   =   meshgrid(lon_res2,lat_res2,depth_res2);
mVmodv_o                    =    griddata(modvOrCartesian(:,1),...
    modvOrCartesian(:,2),modvOrCartesian(:,3),modvOrCartesian(:,4),...
    Xmodv_o,Ymodv_o,Zmodv_o);

sliceY                      =   mean(lat_res2);
sliceX                      =   mean(lon_res2);
sliceZ                      =   mean(depth_res2);

subplot(2,2,2)
slice(Xmodv_o,Ymodv_o,Zmodv_o, mVmodv_o,sliceX,sliceY,sliceZ)
colormap(color)
view(0,90)
hcb                         =   colorbar;
title(hcb,'V (km/s)')
title('Cartesian velocity model')

%% now show inversion velocity model

subplot(2,2,3)
slice(XEqSpace,YEqSpace,ZEqSpace,modvEqS,sliceX,sliceY,sliceZ)
colormap(color)
view(0,90)
hcb                         =   colorbar;
title(hcb,'V (km/s)')
title('Equally-spaced velocity model')

%% propagation model

subplot(2,2,4)
set(gca,'zdir','rev')
slice(Xq,Yq,Zq,modvP,sliceY,sliceX,sliceZ)
colormap(color)
view(0,90)
hcb                         =   colorbar;
title(hcb,'V (km/s)')
title('Propagation velocity model')

storeFolder                     =   'Tests';
FName_testV                     =   'VTest';
pathFolder                      =...
    fullfile(FPath,FLabel,storeFolder,FName_testV);
saveas(testModv, pathFolder);    

end