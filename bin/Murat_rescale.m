function [X,Y,Z,modvR,mR]	=   Murat_rescale(x_o,y_o,z_o,v_o,x,y,z)
% function [modvR,mR]	=   Murat_rescale(x_o,y_o,z_o,v_o,x,y,z)
% 
% RESCALE field to output reference system
%
% Input parameters:
%    x_o,y_o,z_o:	original reference system
%    v_o:           original field
%
% Output parameters:
%    modvR:         rescaled field

[X,Y,Z]                     =   meshgrid(x,y,z);
mR                          =   griddata(x_o,y_o,z_o,v_o,X,Y,Z);

if find(isnan(mR))
    mR             =   inpaintn(mR);
end

modvR                       =   Murat_unfold(X,Y,Z,mR);
end

    