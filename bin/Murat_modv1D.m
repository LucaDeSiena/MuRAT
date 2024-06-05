function [modv,pvel,modvPlot]...
                            =...
                            Murat_modv1D(modvXYZ,modvOriginal,PorS,origin)
% function [modv,pvel,modvPlot]=...
%                           Murat_modv1D(modvXYZ,modvOriginal,PorS,origin)
% CREATES the velocity models when using a 1D input
%
% Input parameters:         
%    modvXYZ:               output grids of MuRAT_fold
%    modvOriginal:          original 1D velocity model
%    PorS:                  choose if P- or S-waves for ray measurements
%    origin:                origin of the model
%
% Output parameters:
%    modv:                  propagation and inversion velocity model
%    pvel:                  propagation velocity model in matrix
%    modvPlot:              velocity model to be plot

modv                        =   modvXYZ;
z1D                         =   -modvOriginal(:,1)*1000;
v1D                         =   modvOriginal(:,PorS+1);

for k=1:length(modv(:,1))
    zModv_k                 =   modv(k,3);
    [~,indexModv]           =   min(abs(z1D-zModv_k));
    modv(k,4)               =   v1D(indexModv);
end
    
% Nodes of the original velocity model
xD                          =   unique(modv(:,1));
yD                          =   unique(modv(:,2));
zD                          =   sort(unique(modv(:,3)),'descend');

% Create its meshgrid
[~,~,~,pvel]                =   Murat_fold(xD,yD,zD,modv(:,4));

% Change to lat/lon for plotting
wgs84                       =   wgs84Ellipsoid("m");
d                           =   sqrt(modv(:,1).^2 + modv(:,2).^2);
az                          =   atan(modv(:,1)./modv(:,2))*360/2/pi;
az(isnan(az))               =   0;
[lat2,lon2]                 =   reckon(origin(1),origin(2),d,az,wgs84);
modvPlot                    =   [lon2 lat2 modv(:,3:4)];

end
