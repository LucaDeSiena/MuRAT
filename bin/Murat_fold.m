function [X,Y,Z,VMesh]          =   Murat_fold(x,y,z,v)
% function [X,Y,Z,VMesh]          =   Murat_fold(x,y,z,v)
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

lx                          =   length(x);
ly                          =   length(y);
lz                          =   length(z);
VMesh                       =   zeros(ly,lx,lz);

[X,Y,Z]                     =   meshgrid(x,y,z);
    
if (nargin==4)
    index                   =   0;
    for i=1:lx
        for j=1:ly
            for k=1:lz
                index       =   index+1;
                VMesh(j,i,k)=   v(index);
            end
        end
    end
end

end

    