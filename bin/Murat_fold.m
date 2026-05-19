function [X,Y,Z,VMesh]      =   Murat_fold(x,y,z,v)
% function [X,Y,Z,VMesh]    =   Murat_fold(x,y,z,v)
%
% MESHGRID 3D for many fields used by MuRAT, including meshgrid switch
%
% Input parameters:
%    x:         x vector
%    y:         y vector
%    z:         z vector
%    v:         field vector - optional
%
% Output parameters:
%    X:         3D x matrix in meshgrid format
%    Y:         3D y matrix in meshgrid format
%    Z:         3D z matrix in meshgrid format
%    VMesh:     3D field matrix in meshgrid format

x       =   x(:).';
y       =   y(:).';
z       =   z(:).';

lx      =   numel(x);
ly      =   numel(y);
lz      =   numel(z);

[X,Y,Z] =   meshgrid(x,y,z);

% Ensure the field vector v is provided; if not, initialize VMesh to empty
if nargin < 4
    VMesh   =   [];
else
    VMesh   =   permute(reshape(v, [lz, ly, lx]), [2,3,1]);
end
end

    