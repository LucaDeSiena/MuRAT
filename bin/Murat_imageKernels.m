function Murat_imageKernels(X,Y,Z,V,color,sections)
% function Murat_imageKernels(X,Y,Z,V,color,sections)
%
% PLOTS a 3D image for the kernels
%
% Input parameters:
%    X:         3D x matrix
%    Y:         3D y matrix
%    Z:         3D z matrix
%    V:         3D field matrix
%    color:     name of the colormap
%    sections:  location of sections   

V(isinf(V))             =   10^-100;
slice(X, Y, Z, V, sections(1), sections(2), sections(3))
shading flat; colormap(color); colorbar
grid on
ax                      =   gca;
ax.GridLineStyle        =   '-';
ax.GridColor            =   'k';
ax.GridAlpha            =   1;
ax.LineWidth            =   1;

