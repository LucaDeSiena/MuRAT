function Murat                  =   Murat_checks(Murat)
%ADDITIONAL input variables that are not set by the user.

disp('Checks')

% Creating folder to store results
if exist(cat(2,'./',Murat.input.label),'dir')~=7
    
    mkdir(cat(2,'./',Murat.input.label))
    mkdir(cat(2,'./',...
        Murat.input.label,'/Rays_Checks'))
    mkdir(cat(2,'./',...
        Murat.input.label,'/Results'))
    mkdir(cat(2,'./',...
        Murat.input.label,'/Resolution'))
    mkdir(cat(2,'./',...
        Murat.input.label,'/VTK'))
    mkdir(cat(2,'./',...
        Murat.input.label,'/TXT'))
    
end

% Get general paths/data options
[Murat.input.listSac,~]         =...
    createsList(Murat.input.dataDirectory);

if Murat.input.nonLinear == 1
    
    % Number of time windows to compute coda intensity
    Murat.input.fitNumber       =...
        Murat.input.codaWindow/Murat.input.fitLengthWindows;
    
    Murat.input.fitTrialQc      =...
        (Murat.input.minimumQcGrid +...
        (Murat.input.maximumQcGrid - Murat.input.minimumQcGrid)/...
        (Murat.input.totalQcNumber-1)*(0:Murat.input.totalQcNumber-1))';
    
end

%% VELOCITY MODELS: ORIGINAL, INVERSION, and PROPAGATION
% Loading inputs to create the additional parameters
origin                          =   Murat.input.origin;
ending                          =   Murat.input.end;
nLat                            =   Murat.input.gridLat;
nLong                           =   Murat.input.gridLong;
nzc                             =   Murat.input.gridZ;
velocityModel                   =   ['velocity_models/',Murat.input.namev];

% Save x,y,z in degrees switching as longitude comes second
Murat.input.x                   =   linspace(origin(2),ending(2),nLong);
Murat.input.y                   =   linspace(origin(1),ending(1),nLat);
Murat.input.z                   =   linspace(origin(3),ending(3),nzc);

% Find distance and azimuth to change in meters - requires longitude first
dist_xdeg                       =   ending(2)-origin(2);
dist_ydeg                       =   ending(1)-origin(1);
dist_x                          =   deg2km(dist_xdeg)*1000;
dist_y                          =   deg2km(dist_ydeg)*1000;

% Coordinates of inversion points in meters
xMeters                         =   linspace(0,dist_x,nLong);
yMeters                         =   linspace(0,dist_y,nLat);
zMeters                         =   linspace(origin(3),ending(3),nzc);

Murat.input.gridStepX           =   xMeters(2)-xMeters(1);
Murat.input.gridStepY           =   yMeters(2)-yMeters(1);

if Murat.input.availableVelocity ==  0
    % The velocity model is a false 3D, it is used as
    % both inversion and propagation model
    modv                        =...
        Murat_unfold(xMeters',yMeters',zMeters');
    modv1D                      =   load(velocityModel);
    zIasp91                     =   -modv1D(:,1)*1000;
    
    if Murat.input.POrS == 2
        vIasp91                 =   modv1D(:,3);
    elseif Murat.input.POrS == 3
        vIasp91                 =   modv1D(:,4);
    end
    
    for k=1:length(modv(:,1))
        zModv_k                 =   modv(k,3);
        [~,indexModv]           =   min(abs(zIasp91-zModv_k));
        modv(k,4)               =   vIasp91(indexModv);
    end
    
    % Grid for ray tracing
    gridPropagation.x           =   xMeters';
    gridPropagation.y           =   yMeters';
    gridPropagation.z           =   zMeters';
    Murat.input.gridPropagation =   gridPropagation;
    
    % Nodes of the original velocity model
    xD                          =   unique(modv(:,1));
    yD                          =   unique(modv(:,2));
    zD                          =   sort(unique(modv(:,3)),'descend');
    
    % Create its meshgrid
    [~,~,~,pvel]                =   Murat_fold(xD,yD,zD,modv(:,4));
    Murat.input.modv            =   modv;
    Murat.input.modvp           =   modv;
    Murat.input.modvPlot        =   [];
    
elseif Murat.input.availableVelocity ==  1
    
    % This works in [x,y,z][m]. If you have a 1D you need
    % to prepare a 3D from it spanning the horizontal space,
    
    % Original 3D velocity model from text file.
    modv_o                      =   load(velocityModel);
    
    % Nodes of the original velocity model
    xD                          =   unique(modv_o(:,1));
    yD                          =   unique(modv_o(:,2));
    zD                          =   sort(unique(modv_o(:,3)),'descend');
    
    % Create its meshgrid
    [XD,YD,ZD,V]                =   Murat_fold(xD,yD,zD,modv_o(:,4));
    
    % Interpolated axes for inversion model
    [X,Y,Z]                     =   meshgrid(yMeters,xMeters,zMeters);
    
    % Interpolated 3D inversion model
    mV                          =   interp3(XD,YD,ZD,V,X,Y,Z);
    
    %In case limits outside of the grid interpolate better
    if find(isnan(mV))
        mV                      =   inpaintn(mV);
    end
    Murat.input.modvPlot        =   mV;
    
    %Create the classical velocity model for inversion
    Murat.input.modv             =...
        Murat_unfold(xMeters',yMeters',zMeters',mV);
    
    %Interpolated model for ray-tracing - half grid step for original
    resol2x                     =   abs(xD(2)-xD(1))/2;
    resol2y                     =   abs(yD(2)-yD(1))/2;
    resol2z                     =   abs(zD(2)-zD(1))/2;
    
    % Interpolated vectors
    xp                          =   xD(1):resol2x:xD(end);
    yp                          =   yD(1):resol2y:yD(end);
    zp                          =   zD(1):-resol2z:zD(end);
    
    % Grid for ray tracing
    gridPropagation.x           =   xp;
    gridPropagation.y           =   yp;
    gridPropagation.z           =   zp;
    Murat.input.gridPropagation =   gridPropagation;
    
    % Meshes for interpolation
    [Xq,Yq,Zq]                  =   meshgrid(yp,xp,zp);
    
    % pvel is the propagation velocity 3D grid
    pvel                        =   interp3(XD,YD,ZD,V,Xq,Yq,Zq);
    
    % Standard unpacking of the velocity model for propagation
    modvp                       =   Murat_unfold(xp',yp',zp',pvel);
    Murat.input.modvp           =   modvp(:,1:4);
    
end

Murat.input.pvel                =   pvel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% DEVELOPMENT: PEAK DELAY INVERSION
% % Input kappa value
% vK                            = Murat.input.vonKarman;
%
% % Sets the three parameters for the modelling of peak delay (Takahashi et
% % al. 2008, Table 1).
% vK_x                          = 0.1:0.1:1;
%
% bp_x                          =...
%     [0.1;0.17;0.23;0.28;0.31;0.34;0.36;0.37;0.37;0.38];
%
% C_x                           =...
%     [0.56;1.06;1.56;2;2.28;2.31;2.14;1.9;1.68;1.5];
%
% p_x                             =...
%     [1.19;1.38;1.56;1.71;1.83;1.91;1.95;1.98;1.99;1.99];
%
% F_bp_x                          = griddedInterpolant(vK_x,bp_x);
%
% F_C_x                           = griddedInterpolant(vK_x,C_x);
%
% F_p_x                           = griddedInterpolant(vK_x,p_x);
%
% Murat.input.Bessel               = [F_bp_x(vK) F_C_x(vK) F_p_x(vK)];
%

function [listWithFolder,listNoFolder]  =   createsList(directory)
%CREATES a list of visible files in a folder, outputs both with and
% without folder

% Get general paths/figure options
list                                    =   dir(directory);
list                                    =...
    list(~startsWith({list.name}, '.'));

listWithFolder                          =...
    fullfile({list.folder},{list.name})';

listNoFolder                            =   {list.name}';
