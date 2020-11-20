function Ac_i                   =...
    Murat_codaMatrix(Murat,flag,tCoda_i,eventStation_i)
%CREATES the coda attenuation inversion matrix

tWm                             =   Murat.input.codaWindow;
kT                              =   Murat.input.kernelTreshold;
modv                            =   Murat.input.modv;
vS                              =   Murat.input.averageVelocityS;
origin                          =   Murat.input.origin;
sections                        =   Murat.input.sections;

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

[X,Y,Z,~]                       =   Murat_fold(x,y,z);

% Kernel in inversion grid space
mK                              =   interp3(Xk,Yk,Zk,K,X,Y,Z);

%In case limits outside of the grid interpolate better
if find(isnan(mK))
    mK                          =   inpaintn(mK);
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
