function [modv,pvel]        =   Murat_modv1D(modvXYZ,modvOriginal,PorS)
% function [modv,pvel]      =  Murat_modv1D(modvXYZ,modvOriginal,PorS)
% CREATES the velocity models when using a 1D as input
%
% Input parameters:         
%    modvXYZ:               output grids of MuRAT_fold
%    modvOriginal:          original 1D velocity model
%    PorS:                  choose if P- or S-waves for ray measurements
%
% Output parameters:
%    modv:                  propagation and inversion velocity model
%    pvel:                  propagation velocity model in matrix

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
end