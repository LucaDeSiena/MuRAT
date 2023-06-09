% ADDITIONAL input variables that are not set by the user.
function Murat                  =   Murat_checks(Murat)

% INPUTS
dataDirectory                   =   ['./' Murat.input.dataDirectory];
FPath                           =   './';
FLabel                          =   Murat.input.label;
PTime                           =   ['SAChdr.times.' Murat.input.PTime];
PorS                            =   Murat.input.POrS;

origin                          =   Murat.input.origin;
ending                          =   Murat.input.end;
nLat                            =   Murat.input.gridLat;
nLong                           =   Murat.input.gridLong;
nzc                             =   Murat.input.gridZ;
availableVelocity               =   Murat.input.availableVelocity;
velocityModel                   =   ['velocity_models/',Murat.input.namev];

if isempty(Murat.input.originTime)
    originTime                  =   'SAChdr.times.o';
else
    originTime                  =   ['SAChdr.times.' Murat.input.originTime];
end

if isempty(Murat.input.STime)
    STime                       =   'SAChdr.times.t0';
else
    STime                       =   ['SAChdr.times.' Murat.input.STime];
end

Murat.input.originTime          =   originTime;
Murat.input.PTime               =   PTime;
Murat.input.STime               =   STime;


if exist('./temp','dir')==7    
    delete('./temp/*')
else
    mkdir('./temp')
end

%% Creating folders to store results
if exist(FLabel,'dir')==7
    rmdir(FLabel,'s')
end

mkdir(FLabel)
mkdir([FLabel,'/RaysKernels'])
mkdir([FLabel,'/Tests'])
mkdir([FLabel,'/Results'])
mkdir([FLabel,'/Results/PeakDelay'])
mkdir([FLabel,'/Results/Qc'])
mkdir([FLabel,'/Results/Q'])
mkdir([FLabel,'/Results/Parameter'])
mkdir([FLabel,'/Checkerboard'])
mkdir([FLabel,'/Checkerboard/Qc'])
mkdir([FLabel,'/Checkerboard/Q'])
mkdir([FLabel,'/Spike'])
mkdir([FLabel,'/Spike/Qc'])
mkdir([FLabel,'/Spike/Q'])
mkdir([FLabel,'/TXT'])

%% Checking data
[Murat.input.listSac,~]         =   createsList([dataDirectory '/*.sac']);
[header,flag]                   =...
    Murat_testData(dataDirectory,originTime,PTime,STime);
mLat                            =...
    [min(cell2mat(header{:,5})) max(cell2mat(header{:,5}))];
mLon                            =...
    [min(cell2mat(header{:,6})) max(cell2mat(header{:,6}))];

if isequal(flag,1)
    warning('Missing origin times.')
end

if isequal(flag,2)
    warning('Missing S-wave times.')
end
%% VELOCITY MODELS: ORIGINAL, INVERSION, and PROPAGATION
% Save x,y,z in degrees switching as longitude comes second
% Find distance and azimuth to change in meters

modvOriginal                    =   load(velocityModel);

wgs84                           =   wgs84Ellipsoid("m");
[d,az]                          =...
    distance(origin(1),origin(2),ending(1),ending(2),wgs84);
dist_x                          =   d*sin(az*2*pi/360);
dist_y                          =   d*cos(az*2*pi/360);

% Coordinates of inversion points in meters
xM                              =   linspace(0,dist_x,nLong)';
yM                              =   linspace(0,dist_y,nLat)';
zM                              =   linspace(origin(3),ending(3),nzc)';
modvXYZ                         =   Murat_unfoldXYZ(xM,yM,zM);

Murat.input.x                   =   linspace(origin(2),ending(2),nLong)';
Murat.input.y                   =   linspace(origin(1),ending(1),nLat)';
Murat.input.z                   =   linspace(origin(3),ending(3),nzc)';
Murat.input.gridStepX           =   xM(2)-xM(1);
Murat.input.gridStepY           =   yM(2)-yM(1);

if availableVelocity ==  0

    gridPropagation.x           =   xM';
    gridPropagation.y           =   yM';
    gridPropagation.z           =   zM';
    
    [modv,pvel,modvPlot]        =...
        Murat_modv1D(modvXYZ,modvOriginal,PorS,origin);
    
    Murat.input.modv            =   modv;
    Murat.input.modvp           =   modv;
    Murat.input.modvPlot        =   modvPlot;
    
    
elseif availableVelocity ==  1
    
    [modvP,pvel,modvEqS,modvPlo]=   Murat_modv3D(FPath,FLabel,...
        modvOriginal,origin,mLat,mLon,nLat,nLong,nzc);
    
    gridPropagation.x           =   unique(modvP(:,1));
    gridPropagation.y           =   unique(modvP(:,2));
    gridPropagation.z           =   sort(unique(modvP(:,3)),'descend');
    
    Murat.input.modv            =   modvEqS;
    Murat.input.modvp           =   modvP;    
    Murat.input.modvPlot        =   modvPlo;
    
end

Murat.input.gridPropagation     =   gridPropagation;
Murat.input.pvel                =   pvel;
Murat.input.header              =   header;
DD_coord                        =   Murat_DDcoordinates(origin,ending,nLat,nLong,nzc);
Murat.input.DDcoordinates       =   DD_coord;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [listWithFolder,listNoFolder]...
                                =   createsList(directory)
% CREATES a list of visible files in a folder, outputs both with and
% without folder

list                            =   dir(directory);
list                            =   list(~startsWith({list.name}, '.'));

listWithFolder              =	fullfile({list.folder},{list.name})';
listNoFolder                =   {list.name}';
