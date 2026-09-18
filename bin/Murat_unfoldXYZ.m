function r      =   Murat_unfoldXYZ(x,y,z)
% function r    =   Murat_unfold(x,y,z)
%
% ACCEPTS vertical vectors and unfolds them in standard format
%
% Input parameters:
%    x:         x vector
%    y:         y vector
%    z:         z vector
%
% Output parameters:
%    r:         3D field in Murat format

x       =   x(:);
y       =   y(:);
z       =   z(:);

lx      =   numel(x);
ly      =   numel(y);
lz      =   numel(z);

% x: repeat each x for ly*lz entries
xcol    =   repelem(x, ly * lz);

% y: repeat each y for lz entries, then tile for lx blocks
yblock  =   repelem(y, lz);
ycol    =   repmat(yblock, lx, 1);

% z: tile z for each x-y pair
zcol    =   repmat(z, lx * ly, 1);

r       =   [xcol, ycol, zcol];

end