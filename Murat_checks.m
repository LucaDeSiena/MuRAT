function Murat                  =   Murat_checks(Murat)
%ADDITIONAL input variables that are not set by user.

% Creating folder for storing results
if exist(cat(2,Murat.input.workingDirectory,Murat.input.label),'dir')~=7
    mkdir(cat(2,Murat.input.workingDirectory,Murat.input.label))
end

% Check that user has the Mapping Toolbox installed.
Murat.input.hasMappingToolbox   =   license('test', 'map_toolbox');

if ~Murat.input.hasMappingToolbox
    % User does not have the toolbox installed or activated.
    sprintf('No Mapping Toolbox - Maps will not be geo-localised');
end

% Get general paths/data options
[Murat.input.listSac,lista]    =...
    Murat_createsList(Murat.input.dataDirectory); 

if Murat.input.nonLinear == 1
    
    % Number of time windows to compute coda intensity
    Murat.input.fitNumber       =...
        Murat.input.codaWindow/Murat.inversion.fitL;

    Murat.input.fitTrialQc      =...
        (m1min+(m1max-m1min)/(total-1)*(0:total-1))';

end

% DEPRECATED BUT NECESSARY FOR NOW - importing locations from files
if Murat.input.importLocation == 1
    
    Murat.input.Locations       =...
        importLocation('eventLocations.txt','stationLocations.txt',lista);
    
end

% VELOCITY MODELS: ORIGINAL, INVERSION, and PROPAGATION
% Loading inputs to create the additional parameters
origin                          =   Murat.input.origin;
ending                          =   Murat.input.end;
nxc                             =   Murat.input.gridX;
nyc                             =   Murat.input.gridY;
nzc                             =   Murat.input.gridZ;

% Save x,y,z in degrees switching as long comes second
Murat.input.x                   =   linspace(origin(2),ending(2),nxc);
Murat.input.y                   =   linspace(origin(1),ending(1),nyc);
Murat.input.z                   =   linspace(origin(3),ending(3),nzc);

% Find distance and azimuth to change in meters - requires longitude first
dist_xdeg                       =   ending(2)-origin(2);
dist_ydeg                       =   ending(1)-origin(1);
dist_x                          =   deg2km(dist_xdeg)*1000;
dist_y                          =   deg2km(dist_ydeg)*1000;

% Coordinates of inversion points in meters
xMeters                         =   linspace(0,dist_x,nxc);
yMeters                         =   linspace(0,dist_y,nyc);
zMeters                         =   linspace(origin(3),ending(3),nzc);

Murat.input.gridStepX           =   xMeters(2)-xMeters(1);
Murat.input.gridStepY           =   yMeters(2)-yMeters(1);

if Murat.input.availableVelocity ==  0
    % The "velocity model" is simply a 3D with V0 everywhere
    if Matlab.input.PorS == 2
        v                       =   vP;
    elseif Matlab.input.PorS == 3
        v                       =   vS;
    end
    % Inversion and propagation model
    modv                        =...
        Murat_unfold(xMeters',yMeters',zMeters');
    modv(:,4)                   =   v;
    
    Murat.input.modv            =   modv;
    Murat.input.modvp           =   modv;
    
    % Grid for ray tracing
    Murat.input.gridD           =...
        createsGrid(xMeters,yMeters,zMeters);
    
elseif Murat.input.availableVelocity ==  1
        
    % This works in [x,y,z][m]. If you have a 1D you need
    % to prepare a 3D from it spanning the horizontal space,
    
    % Original 3D velocity model from text file.
    modv_o                      =   load(Murat.input.namev);
    modv_o(:,5)                 =   0;
    
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
    
    % Meshes for interpolation
    [Xq,Yq,Zq]                  =   meshgrid(yp,xp,zp);
    
    % pvel is the propagation velocity 3D grid
    pvel                        =   interp3(XD,YD,ZD,V,Xq,Yq,Zq);
    Murat.input.pvel            =   pvel;
    
    % Standard unpacking of the velocity model for propagation
    modvp                       =   Murat_unfold(xp',yp',zp',pvel);
    modvp(:,5)                  =   0;
    Murat.input.modvp           =   modvp;
    
    % Grid for ray tracing
    Murat.input.gridD           =    createsGrid(xp,yp,zp);
    
end

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

function evestaz                =   importLocation(evenFile,stazFile,list)
%LOADS event and station coordinates from two files. DEPRECATED as it
%requires mathing the names of the external files.

evestaz                         =   zeros(length(list),6);

%Name of the event file if importing event locations from file
namee                           =   evenFile;

% Opening the event file. Format is:
% column (1) = twelve numbers for the origin time of the event
% (date+time in seconds)
% column (2) = latitude
% column (3) = longitude
% column (4) = Altitude above sea level in meters
event                           =   fopen(namee);
namee                           =   textscan(event,'%s %f %f %f');
nameeven                        =   namee{1};
even                            =   [namee{2} namee{3} -namee{4}];
fclose(event);

%Name of the station file if importing event locations from file
names                           =   stazFile;

% Opening the station file. Format is:
% column (1) = Name of station (3 characters)
% column (2) = latitude
% column (3) = longitude
% column (4) = Altitude above sea level in meters
station                         =   fopen(names);
names                           =   textscan(station,'%s %f %f %f');
namestation                     =   names{1};
stat                            =   [names{2} names{3} names{4}];
fclose(station);

indexray=0;
for i = 1:length(list)
    %Here it takes the event/station name. It is necessary to adapt the
    %numbers to where event name and station name are.
    li                          =   list{i};
    li1                         =   cat(2,li(1:12),li(18:20));
    
    for ii = 1:length(nameeven)
        namee1                  =   nameeven{ii};
        namee                   =   namee1(1:12);
        
        % Loop over stations in staz.txt file
        for ir = 1:length(namestation)
            namest1             =   namestation{ir};
            namest              =   namest1(1:3);
            levst               =   cat(2,namee,namest);
            
            if find(strncmp(li1,levst,length(levst)))>0
                indexray        =   indexray+1;
                sst             =   [even(ii,:) stat(ir,:)];
                evestaz(indexray,1:6) = sst;
                break
            end
        end
    end
    
end

function gridD                  =   createsGrid(x,y,z)

lx                              =   length(x);
ly                              =   length(y);
lz                              =   length(z);

l                               =   max([lx ly lz]);

% Same
gridD                           =   zeros(3,max(l));
gridD(1,1:lx)                   =   x;
gridD(2,1:ly)                   =   y;
gridD(3,1:lz)                   =   z;
    
    
