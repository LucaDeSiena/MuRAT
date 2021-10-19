function modv               =   Murat_unfold(X,Y,Z,mV)
% function modv                  =   Murat_unfold(X,Y,Z,mV)
%
% ACCEPTS matrices created by meshgrid and unfolds them in MuRAT format
%
% Input parameters:
%    x:         X matrix from meshgrid
%    y:         Y matrix from meshgrid
%    z:         Z matrix from meshgrid
%    v:         V matrix from meshgrid
%
% Output parameters:
%    modv:      field in Murat format

ly =size(X,1);
lx =size(X,2);
lz =size(X,3);
modv=zeros(lx*ly*lz,4);
if (nargin==4)
    index                   =   0;
    for i=1:lx
        for j=1:ly
            for k=1:lz
                index       =   index+1;                
                modv(index,1:4)...
                            =   [X(j,i,k) Y(j,i,k) Z(j,i,k) mV(j,i,k)];
            end
        end
    end
end

% Line to get it in Murat format;
end
