function AQc_i    =   Murat_codaMatrix(modvQc,K_grid_start,r_grid_start,...
    K_grid_end,r_grid_end,flag,origin,sections,FName,pathFolder)
% function AQc_i  =   Murat_codaMatrix(modvQc,K_grid_start,r_grid_start,...
%   K_grid_end,r_grid_end,flag,origin,sections)
%
% CREATES the coda attenuation inversion matrix and plots the corresponding
%   kernels
%
% Input parameters:         
%    modvQc:                velocity model for grid
%    K_grid:                kernel from Murat_kernels
%    r_grid:                grid from Murat_kernels
%    flag:                  if 1 creates the figure for first couple
%    origin:                origin of the grid
%    sections:              sections of the figure
%
% Output parameters:
%    AQc_i:                 coda inversion matrix

% Nodes of the kernel model space
xKStart             =   unique(r_grid_start(:,1));
yKStart             =   unique(r_grid_start(:,2));
zKStart             =   sort(unique(r_grid_start(:,3)),'descend');
[XkS,YkS,ZkS,KStart]=   Murat_fold(xKStart,yKStart,zKStart,K_grid_start);

xKEnd               =   unique(r_grid_end(:,1));
yKEnd               =   unique(r_grid_end(:,2));
zKEnd               =   sort(unique(r_grid_end(:,3)),'descend');
[XkE,YkE,ZkE,KEnd]  =   Murat_fold(xKEnd,yKEnd,zKEnd,K_grid_end);

% Interpolated axes for inversion model
x                   =   unique(modvQc(:,1));
y                   =   unique(modvQc(:,2));
z                   =   sort(unique(modvQc(:,3)),'descend');
[X,Y,Z,~]           =   Murat_fold(x,y,z);

% Kernels in inversion grid space
mKS                 =   interp3(XkS,YkS,ZkS,KStart,X,Y,Z);
mKE                 =   interp3(XkE,YkE,ZkE,KEnd,X,Y,Z);

% Replace NaNs and zeros with a tiny positive floor; if entire grid empty,
% create a single-point kernel at the nearest grid node to the max original
tiny = realmin;            % machine smallest positive
fixKernel           =   @(m, K_grid, r_grid, xg, yg, zg) ...
    fixKernelLocal(m, K_grid, r_grid, xg, yg, zg, tiny);

if any(isnan(mKS(:))) || all(mKS(:) == 0)
    mKS             =   fixKernel(mKS,K_grid_start,r_grid_start,x,y,z);
end
if any(isnan(mKE(:))) || all(mKE(:) == 0)
    mKE             =   fixKernel(mKE,K_grid_end,r_grid_end,x,y,z);
end

mK                  =   mKE-mKS;
mK                  =   fixKernel(mK,K_grid_end,r_grid_end,x,y,z);

% Kernel in its grid space
if flag == 1
    % The next figure checks the sensitivity of coda attenuation
    % measurements. The code creates figures that show sections in the
    % sensitivity kernels. The left panel shows the sensitivity kernel 
    % in the full space while the rigth panel shows the normalized
    % kernel in the inversion grid.
    
    % reorder sections and compute geographic sampling
    sections1       =   [sections(2) sections(1) sections(3)];
    d               =   sqrt(X.^2 + Y.^2);
    az              =   atan2d(X, Y);
    az(isnan(az))   =   0;
    wgs84           =   wgs84Ellipsoid("m");
    [Y1,X1]         =   reckon(origin(1),origin(2),d,az,wgs84);

    % query grid in geographic coords and depth
    x1              =   linspace(min(X1(:)),max(X1(:)),length(x))';
    y1              =   linspace(min(Y1(:)),max(Y1(:)),length(y))';
    z1              =   z;
    [XEqS,YEqS,ZEqS]=   meshgrid(x1,y1,z1);

    % query grid in geographic coords and depth
    mEqSpace        =   griddata(X1,Y1,Z,mK,XEqS,YEqS,ZEqS);

    % kernel grid to geographic coords
    Xk1             =   origin(2) + km2deg(XkS/1000);
    Yk1             =   origin(1) + km2deg(YkS/1000);
    
    kernels         =   myfig(FName);
    
    subplot(1,2,1)
    Murat_imageKernels(Xk1,Yk1,ZkS/1000,KStart,inferno,sections1)
    
    subplot(1,2,2)
    Murat_imageKernels(XEqS,YEqS,ZEqS,mEqSpace,inferno,sections1)

    Murat_saveFigures_2panels(kernels,pathFolder);
end

%pre-define 3D matrix in space
AQc_i               =   reshape(permute(mK, [3 1 2]), [], 1);
% posMask             =   AQc_i > 0;
% if any(posMask)
%     minPos          =   min(AQc_i(posMask));
%     AQc_i(AQc_i < 0)=   minPos;
% end

% Residual from cutting the grid (it is always < 1%).
AQc_i               =   AQc_i/sum(AQc_i);
end

% --- Local helper to repair kernel when interp produced NaNs or zeros ---
function m = fixKernelLocal(m, K_grid, r_grid, xg, yg, zg, tiny)
    % set NaNs and zeros to tiny
    m(isnan(m))     =   tiny;
    m(m == 0)       =   tiny;

    % if still effectively empty (all tiny), place a single 1 at nearest node
    if ~any(m > tiny)
        modGrid     =   Murat_unfoldXYZ(xg, yg, zg);   % N x 3
        % pick the original kernel's max location
        [~, idxMax] =   max(K_grid);
        rmax        =   r_grid(idxMax, :);
        dists       =   sum((modGrid(:,1:3) - rmax).^2, 2);
        [~, imin]   =   min(dists);
        tmp         =   zeros(size(modGrid,1),1);
        tmp(imin)   =   1;
        [~,~,~,m]   =   Murat_fold(xg, yg, zg, tmp);
    end
end
