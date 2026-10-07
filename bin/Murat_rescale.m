function [X, Y, Z, modvR, mR] = Murat_rescale(x_o, y_o, z_o, v_o, x, y, z)
% MURAT_RESCALE  Interpolate a field onto a new regular grid.
%
%   [X,Y,Z,modvR,mR] = Murat_rescale(x_o,y_o,z_o,v_o,x,y,z)
%
%   Input parameters:
%       x_o, y_o, z_o   original (scattered or regular) coordinate vectors
%       v_o             original field values at those coordinates
%       x, y, z         target grid coordinate vectors
%
%   Output parameters:
%       X, Y, Z         meshgrid arrays of the target grid
%       mR              interpolated field on the target grid (no NaNs)
%       modvR           mR unfolded to a column vector via Murat_unfold
%
%   griddata() can return NaN for target points that lie outside the convex
%   hull of the source points.  Those gaps are filled by nearest-neighbour
%   propagation along each grid dimension — a base-MATLAB replacement for
%   the GIBBON inpaintn() routine (which carried an uncertain licence).

[X, Y, Z]   =   meshgrid(x, y, z);
mR          =   griddata(x_o, y_o, z_o, v_o, X, Y, Z);

if any(isnan(mR(:)))
    mR      =   fillNaN3(mR);
end

modvR       =   Murat_unfold(X, Y, Z, mR);
end

% -------------------------------------------------------------------------
function mR = fillNaN3(mR)
% FILLNAN3  Fill NaNs in a 3-D array by nearest-neighbour propagation.
%
%   Two passes along each of the three dimensions guarantee that every NaN
%   reachable from a valid value is filled, including the case where an
%   entire border plane is NaN.  Uses only base MATLAB (no toolboxes).

for pass = 1:2
    for dim = 1:3
        mR = fillmissing(mR, 'nearest', dim);
    end
    if ~any(isnan(mR(:))), return; end
end
end
