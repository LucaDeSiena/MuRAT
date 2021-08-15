function [X,Y,Z,V]          =   Murat_fold(x,y,z,v)
% function [X,Y,Z,V]          =   Murat_fold(x,y,z,v)
%
% MESHGRID 3D for many fields used by MuRAT
%
% Input parameters:
%    x:         x vector
%    y:         y vector
%    z:         z vector
%    v:         field vector
%
% Output parameters:
%    X:         3D x matrix
%    Y:         3D y matrix
%    Z:         3D z matrix
%    V:         3D field matrix

lx                          =   length(x);
ly                          =   length(y);
lz                          =   length(z);
V                           =   zeros(lx,ly,lz);

[X,Y,Z]                     =   meshgrid(y,x,z);
    
if (nargin==4)
    index                   =   0;
    for i=1:lx
        for j=1:ly
            for k=1:lz
                index       =   index+1;
                V(i,j,k)    =   v(index);
            end
        end
    end
end
    