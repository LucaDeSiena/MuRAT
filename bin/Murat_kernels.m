function [K_grid,r_grid1]	=...
    Murat_kernels(T,event,station,modv,v,kT,B0,Le_1)
% function [K_grid,r_grid1]	=...
%     Murat_kernels(T,event,station,modv,v,kT,B0,Le_1)
%
% COMPUTATION of kernel integral for an event-station couple in a
%   3D grid (modv format), to be used for coda tomography. There is the
%   option to define the average velocity (v) and divide the original
%   spacing by kT - higher kT brings longer computational time.
%
% Input parameters:
%    T:             lapse time
%    event:         event coordinates
%    station:       station coordinates
%    modv:          velocity model used for setting up the grid
%    v:             average velocity
%    kT:            the computational parameter that speeds up calculations
%    B0:            albedo
%    Le_1:          extinction length
%
% Output parameters:
%    r_grid1:           The coordinates of the grid
%    K_grid:            The kernel estimated at each point of the grid
%
% Structure of the file describing the coordinates:
%     Nx    -> number of x values in the grid
%     Ny    -> number of y values in the grid
%     Nz    -> number of z values in the grid
%     x(1)  -> first x coordinate
%     x(2)  ...
%     x(3)
%     ...
%     x(Nx) -> last x coordinate
%     y(1)  -> first y coordinate
%     y(2)
%     ...
%     y(Ny) -> last y coordinate
%     z(1)  -> first z coordinate
%     z(2)
%     ...
%     z(Nz) -> last z coordinate
%
% Structure of the file with the K function results: 
%     K(1)      1st point for x(1),y(1),z(1)
%     K(2)      2nd point for x(2),y(1),z(1)
%     K(3)      3rd point for x(3),y(1),z(1)
%     ....
%     K(Nx * Ny * Nz)  last point for x(Nx),y(Ny),z(Nz)
% Original creators: Angel De La Torre & Edoardo Del Pezzo
% Adapted for coda attenuation tomography by De Siena.

DT                          =   0.5;
%%
% Initial computations to understand were to compute in space
x1                          =   unique(modv(:,1));
y1                          =   unique(modv(:,2));
z1                          =   sort(unique(modv(:,3)),'descend');
stepgx                      =   x1(2)-x1(1);
stepgy                      =   y1(2)-y1(1);
stepgz                      =   z1(2)-z1(1);

DRx                         =   stepgx/kT;
DRy                         =   stepgy/kT;
DRz                         =   stepgz/kT;
DR                          =   min([DRx DRy abs(DRz)])/1000;

% Origin of the grid: middle point between source and receiver
origin                      =   (event(1:3)+station(1:3))/2;

% Add a cushion of 5 seconds
S                           =   (v*T+5)*1000;

% New reference system
xyzs                        =   event(1:3)-origin(1:3);
xyzr                        =   station(1:3)-origin(1:3);

% Distance between source and receiver
dxyz                        =   xyzr(1:3)-xyzs(1:3);
D0                          =...
    sqrt(dxyz(1)^2 + dxyz(2)^2 + dxyz(3)^2)/1000;

% Grid definition
x_grid                      =   -S+DRx:DRx:S;
y_grid                      =   -S+DRy:DRy:S;
z_grid                      =   max(z1)-origin(3):DRz:-S;

Nx                          =   length(x_grid); 
Ny                          =   length(y_grid); 
Nz                          =   length(z_grid);
Nxyz                        =   Nx*Ny*Nz;

if Nxyz>50e6
    fprintf('Warning: a lot of points are going to be included in the grid (%d points)\n',Nx*Ny*Nz)
    fprintf('Press any key to continue, or CRTL-C to abort\n')
    pause
end

% Sets the 3D matrix
r_grid                      =   Murat_unfold(x_grid',y_grid',z_grid');


% Distances r1 and r2, from each point of the grid to source or to receptor
xyz                         =   r_grid(:,1:3);
d1xyz                       =   xyz-xyzs;
d2xyz                       =   xyz-xyzr;
r1_grid                     =...
    sqrt(d1xyz(:,1).^2 + d1xyz(:,2).^2 + d1xyz(:,3).^2)/1000;
r2_grid                     =...
    sqrt(d2xyz(:,1).^2 + d2xyz(:,2).^2 + d2xyz(:,3).^2)/1000;

% Elements of the grid inside the ellipsoid
cond_interior               =   r1_grid + r2_grid < T*v;

%%
% Catalog of Paasschens Functions

% Vector of distances
R                           =   DR:DR:T*v;

% Max number points coda
NT                          =   ceil((T+DT)/DT); 

% Number of distances
NR                          =   length(R);

% To store t0 for each r
t0_R                        =   zeros(NR,1);

% To store A for each r
A_R                         =   zeros(NR,1);

% Store the number of points in coda for each r
N_R                         =   zeros(NR,1); 

% Store the coda for each r
coda_R                      =   zeros(NR,NT);

for idx_r=1:NR
    r                       =   R(idx_r);
    [t0,A,N,coda]           =...
        Murat_paasschensFunction(r,v,B0,Le_1,DT,T);
    t0_R(idx_r)             =   t0;
    A_R(idx_r)              =   A;
    N_R(idx_r)              =   N;
    coda_R(idx_r,1:N)       =   coda(1:N);
end

%%
% Catalog of r1 - r2 convolutions evaluated at T
R1                          =   R;
R2                          =   R;
coda_coda                   =   zeros(NR,NR);
coda_delta                  =   zeros(NR,NR);

for idx_r1=1:NR
    for idx_r2              =   1:idx_r1
        r1                  =   R1(idx_r1);
        r2                  =   R2(idx_r2);
        
        % delta_T is the overlapping interval between E(r1,t) and E(r2,T-r)
        delta_T             =   T-t0_R(idx_r1) - t0_R(idx_r2);
        
        % Number of samples involved in the overlapping interval
        N_mues              =   round(delta_T/DT);
        
        if r1+r2>=(D0-2*DR) && r1+r2<(T+DT)*v && N_mues>=1
            
            % Computing (coda1 x delta2) and (delta1 x coda 2)
            coda1           =   A_R(idx_r1)*coda_R(idx_r2,N_mues+1); 
            coda2           =   A_R(idx_r2)*coda_R(idx_r1,N_mues+1);
            coda_delta...
                (idx_r1,...
                idx_r2)     =  coda1+coda2;
            
            % Computing the contribution of (coda1 x coda2)
            % integral coda1(u) coda(T-u) du
            x1              =   coda_R(idx_r1,1:N_mues+1);
            x2              =   coda_R(idx_r2,N_mues+1:-1:1);
            coda_coda...
                (idx_r1,...
                idx_r2)     =   sum(x1.*x2)*DT;
            
            % simmetrical calculations
            coda_delta...
                (idx_r2,...
                idx_r1)     =   coda_delta(idx_r1,idx_r2);
            coda_coda...
                (idx_r2,...
                idx_r1)     =   coda_coda(idx_r1,idx_r2);
        end
    end
    
end

K_integral                  =   coda_coda+coda_delta;

% Computing K(x,y,z) by linear interpolation inside the ellipsoid
c2                          =   cond_interior;
K_grid                      =   zeros(Nx*Ny*Nz,1);
K_grid(c2)                  =...
    interp2(R1,R2,K_integral,r1_grid(c2),r2_grid(c2));
K_grid(isnan(K_grid))       =   max(K_grid);
r_grid1                     =...
    [r_grid(:,1)+origin(1) r_grid(:,2)+origin(2) r_grid(:,3)+origin(3)];
