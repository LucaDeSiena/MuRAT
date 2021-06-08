function [X,Y,Z,V]          =   Murat_fold(x,y,z,v)

%3D meshgrid of a field
lx                          =   length(x); % number of positions x
ly                          =   length(y); % number of positions y
lz                          =   length(z); % number of positions z
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
    