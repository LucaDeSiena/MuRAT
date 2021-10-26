% ADDITIONAL input variables that are not set by the user.
function Murat                  =   Murat_checks

[file,path]                     =   uigetfile;

if isequal(file,0)
   error('User selected Cancel!');   
else
   disp(['User selected ', fullfile(path,file)]);
end

run(fullfile(path, file));

% INPUTS
dataDirectory                   =   ['./' Murat.input.dataDirectory]; %#ok<NODEF>
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
    originTime                  =...
        ['SAChdr.times.' Murat.input.originTime];
end

if isempty(Murat.input.STime)
    STime                  =   'SAChdr.times.t0';
else
    STime                  =...
        ['SAChdr.times.' Murat.input.STime];
end

Murat.input.originTime          =   originTime;
Murat.input.PTime               =   PTime;
Murat.input.STime               =   STime;


if exist('./temp','dir')==7    
    delete('./temp/*')
else
    mkdir('./temp')
end

% Checking data
[Murat.input.listSac,~]         =   createsList([dataDirectory '/*.sac']);
[Murat.input.header,flag]       =...
    Murat_testData(dataDirectory,originTime,PTime,STime);

if isequal(flag,1)
    warning('Missing origin times.')
end

if isequal(flag,2)
    warning('Missing S-wave times.')
end
%% VELOCITY MODELS: ORIGINAL, INVERSION, and PROPAGATION
% Save x,y,z in degrees switching as longitude comes second
% Find distance and azimuth to change in meters - requires longitude first
dist_x                          =   deg2km(ending(2)-origin(2))*1000;
dist_y                          =   deg2km(ending(1)-origin(1))*1000;

% Coordinates of inversion points in meters
xM                              =   linspace(0,dist_x,nLong)';
yM                              =   linspace(0,dist_y,nLat)';
zM                              =   linspace(origin(3),ending(3),nzc)';
modvXYZ                         =   Murat_unfoldXYZ(xM,yM,zM);

Murat.input.x                   =   linspace(origin(2),ending(2),nLong);
Murat.input.y                   =   linspace(origin(1),ending(1),nLat);
Murat.input.z                   =   linspace(origin(3),ending(3),nzc);
Murat.input.gridStepX           =   xM(2)-xM(1);
Murat.input.gridStepY           =   yM(2)-yM(1);

modvOriginal                    =   load(velocityModel);
if availableVelocity ==  0

    gridPropagation.x           =   xM';
    gridPropagation.y           =   yM';
    gridPropagation.z           =   zM';
    
    [modv,pvel]                 =...
        Murat_modv1D(modvXYZ,modvOriginal,PorS);
    
    Murat.input.modv            =   modv;
    Murat.input.modvp           =   modv;
    Murat.input.modvPlot        =   [];
    
elseif availableVelocity ==  1
    
    [modvP,modvI,modvIP,pvel]          =...
        Murat_modv3D(modvXYZ,modvOriginal,origin,0);
    
    gridPropagation.x           =   unique(modvP(:,1));
    gridPropagation.y           =   unique(modvP(:,2));
    gridPropagation.z           =   sort(unique(modvP(:,3)),'descend');
    
    Murat.input.modv            =   modvI;
    Murat.input.modvp           =   modvP(:,1:4);    
    Murat.input.modvPlot        =   modvIP;
    
end

Murat.input.gridPropagation     =   gridPropagation;
Murat.input.pvel                =   pvel;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [listWithFolder,listNoFolder]...
                                =   createsList(directory)
% CREATES a list of visible files in a folder, outputs both with and
% without folder

list                            =   dir(directory);
list                            =   list(~startsWith({list.name}, '.'));

listWithFolder                  =	fullfile({list.folder},{list.name})';
listNoFolder                    =   {list.name}';
