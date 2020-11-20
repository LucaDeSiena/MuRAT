function Murat_imageKernels(X,Y,Z,V,color,sections)
%PLOTS a 3D image for the kernels
slice(X, Y, Z, V, sections(1), sections(2), sections(3))
colorbar
shading flat
colormap(color);
colorbar

grid on
    
ax                      =   gca;
ax.GridLineStyle        =   '-';
ax.GridColor            =   'k';
ax.GridAlpha            =   1;
ax.LineWidth            =   1;

