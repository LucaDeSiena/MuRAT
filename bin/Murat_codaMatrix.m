function Ac_i                   =...
    Murat_codaMatrix(Murat,flag,tCoda_i,eventStation_i)
%CREATES the coda attenuation inversion matrix

tWm                             =   Murat.input.codaWindow;
kT                              =   Murat.input.kernelTreshold;
modv                            =   Murat.input.modv;
vS                              =   Murat.input.averageVelocityS;
origin                          =   Murat.input.origin;
sections                        =   Murat.input.sections;

% Coda kernels compute value at the centre of the cell, differently from
% rays, which refer to the shallowest SW corner.
stepgX                          =   (modv(2,1) - modv(1,1))/2;
stepgY                          =   (modv(2,2) - modv(1,2))/2;
stepgZ                          =   (modv(2,3) - modv(1,3))/2;

% Update modv to compute at the centre of the block
modv(:,1)                       =   modv(:,1) + stepgX;
modv(:,2)                       =   modv(:,2) + stepgY;
modv(:,3)                       =   modv(:,3) + stepgZ;

[K_grid,r_grid]                 =...
    Murat_kernels(tCoda_i+tWm/2,eventStation_i(1:3),eventStation_i(4:6),...
    modv,vS,kT);

% Nodes of the kernel model space
xK                              =   unique(r_grid(:,1));
yK                              =   unique(r_grid(:,2));
zK                              =   sort(unique(r_grid(:,3)),'descend');

[Xk,Yk,Zk,K]                    =   Murat_fold(xK,yK,zK,K_grid);

% Interpolated axes for inversion model
x                               =   unique(modv(:,1));
y                               =   unique(modv(:,2));
z                               =   sort(unique(modv(:,3)),'descend');

% Interpolation sets everything in the right place
[X,Y,Z,~]                       =   Murat_fold(x,y,z);

% Kernel in inversion grid space
mK                              =   interp3(Xk,Yk,Zk,K,X,Y,Z);

%In case limits outside of the grid interpolate better
if find(isnan(mK))
    mK(isnan(mK))               =   10^-100;
    mK(mK == 0)                 =   10^-100;
    if isempty(find(mK, 1))
        
        mod_K                   =   Murat_unfold(x,y,z);
        mod_K(:,4)              =   0;
        [~,maxK]                =   max(K_grid);
        rmax                    =   r_grid(maxK,:);
        [~,min_K]               =   min(sqrt((mod_K(:,1)-rmax(1)).^2 ...
        + (mod_K(:,2)-rmax(2)).^2 + (mod_K(:,3)-rmax(3)).^2));
        mod_K(min_K,4)          =   1;
        [~,~,~,mK]              =   Murat_fold(x,y,z,mod_K(:,4));

    end
end

% Kernel in its grid space
if flag == 1
    sections1                   =   [sections(2) sections(1) sections(3)];
    Xk1                         =   origin(2) + km2deg(Xk/1000);
    Yk1                         =   origin(1) + km2deg(Yk/1000);
    X1                          =   origin(2) + km2deg(X/1000);
    Y1                          =   origin(1) + km2deg(Y/1000);

    subplot(1,2,1)
    Murat_imageKernels(Xk1,Yk1,Zk,log(K),'default',sections1)
    
    subplot(1,2,2)
    Murat_imageKernels(X1,Y1,Z,log(mK),'default',sections1)
    
end

%pre-define 3D matrix in space
lx                              =   length(x); % number of positions x
ly                              =   length(y); % number of positions y
lz                              =   length(z); % number of positions z
index                           =   0;
Ac_i                            =   zeros(1,(length(modv(:,1))));
for i=1:lx
    for j=1:ly
        for k=1:lz
            index               =   index+1;
            Ac_i(index)         =   mK(i,j,k);
        end
    end
end

% Residual from cutting the grid (it is always <1%
% resAc = sum(Ac_i)/sum(K_grid);

Ac_i = Ac_i/sum(Ac_i);

%PLOTS a 3D image for the kernels
function Murat_imageKernels(X,Y,Z,V,color,sections)
V(isinf(V))             =   10^-100;

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

