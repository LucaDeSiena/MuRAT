function r                  =   Murat_unfold(x,y,z,mV)
% function r                  =   Murat_unfold(x,y,z,mV)
%
% ACCEPTS vertical vectors and unfolds them in standard format
%
% Input parameters:
%    x:         x vector
%    y:         y vector
%    z:         z vector
%    v:         field vector
%
% Output parameters:
%    r:         field in Murat format

lx                          =   length(x);
ly                          =   length(y);
lz                          =   length(z);
lxyz                        =   lx*ly*lz;

r                           =   zeros(lxyz,3);

index                       =   1;
for i = 1:lx
    index1                  =   index+ly*lz-1;
    gx                      =   repmat(x(i),ly*lz,1);
    r(index:index1,1)       =   gx;
    index                   =   index1+1;
end

index                       =   1;
ry                          =   zeros(ly*lz,1);
for i = 1:ly
    index1                  =   index+lz-1;
    gy                      =   repmat(y(i),lz,1);
    ry(index:index1,1)      =   gy;
    index                   =   index1+1;
end

r(:,2)                      =   repmat(ry,lx,1);
r(:,3)                      =   repmat(z,lx*ly,1);

if (nargin==4)
    index                   =   0;
    for i=1:lx
        for j=1:ly
            for k=1:lz
                index       =   index+1;
                r(index,4)  =   mV(i,j,k);
            end
        end
    end
end
    