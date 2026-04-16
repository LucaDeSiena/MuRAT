% ADDITIONAL checks and input variables that are not set by the user.
function Murat      =   Murat_checks(Murat)
%% Adding paths and inputs
addpath(fullfile(pwd,'bin'));
r       =   fullfile(pwd,'Utilities_Matlab');

% Split genpath into cell array of full folder paths
parts   = strsplit(genpath(r), pathsep);
parts   = parts(~cellfun(@isempty, parts));

% Filters: exclude .git, hidden folders (start with .), and 'private'
isBad = @(p) contains(p, [filesep '.git' filesep], 'IgnoreCase', true) ...
           | contains(p, [filesep '.svn' filesep], 'IgnoreCase', true) ...
           | contains(p, [filesep 'private' filesep], 'IgnoreCase', true) ...
           | any( startsWith( split(p, filesep), '.' ) );

% Keep folders that pass filter and contain at least one .m file
keepMask = false(size(parts));
for k = 1:numel(parts)
    p = parts{k};
    if ~isBad(p)
        % fast check for .m files (not recursive)
        if ~isempty(dir(fullfile(p,'*.m')))
            keepMask(k) = true;
        end
    end
end

addList = parts(keepMask);
if ~isempty(addList)
    addpath(strjoin(addList, pathsep), '-end');
end

if isfolder(Murat.input.label)
    dirOld = Murat.input.label+"_old";
    if isfolder(dirOld)
        rmdir(dirOld,"s");
    else
        movefile(Murat.input.label, dirOld);
    end
    fprintf('Folder "%s" found in current folder. Renamed to folder "%s _old"\n', Murat.input.label,Murat.input.label);
else
    mkdir(Murat.input.label);
end

%% Checking inputs
dataDirectory       =   Murat.input.dataDirectory;
FLabel              =   Murat.input.label;
FPath               =   './';
velocityModel       =   fullfile('velocity_models', Murat.input.namev);

% time fields
PTime               =   ['SAChdr.times.' Murat.input.PTime];
PorS                =   Murat.input.POrS;
origin              =   Murat.input.origin;
ending              =   Murat.input.end;
nLat                =   Murat.input.gridLat;
nLong               =   Murat.input.gridLong;
nzc                 =   Murat.input.gridZ;
availableVelocity   =   Murat.input.availableVelocity;

if isempty(Murat.input.originTime)
    originTime      =   'SAChdr.times.o';
else
    originTime      =   ['SAChdr.times.' Murat.input.originTime];
end

if isempty(Murat.input.STime)
    STime           =   'SAChdr.times.t0';
else
    STime           =   ['SAChdr.times.' Murat.input.STime];
end

Murat.input.originTime  =   originTime;
Murat.input.PTime       =   PTime;
Murat.input.STime       =   STime;

if isfolder('./temp')
    delete(fullfile('./temp','*'))
else
    mkdir('./temp')
end

%% Create folder tree in a loop to avoid repeated mkdir calls
subdirs = { '', 'Rays', 'Kernels', 'Tests', 'Tests/PeakDelay', ...
    'Tests/Qc' 'Tests/Q', 'Tests/LCurve', 'Tests/Clustering', 'Results',...
    'Results/PeakDelay', 'Results/Qc', 'Results/Q','Results/Parameter',...
    'Checkerboard', 'Checkerboard/Qc', 'Checkerboard/PD' ...
    'Checkerboard/Q', 'Spike', 'Spike/Qc', 'Spike/Q', 'TXT' };

if isfolder(FLabel)
    rmdir(FLabel,'s')
end
for i = 1:numel(subdirs)
    mkdir(fullfile(FLabel, subdirs{i}))
end

%% Checking data
[Murat.input.listSac,~] =   createsList([dataDirectory '/*.sac']);
[header,flag,sacHeader] =...
    Murat_testData(dataDirectory,originTime,PTime,STime);

mLat    =   [min(cell2mat(header{:,5})) max(cell2mat(header{:,5}))];
mLon    =   [min(cell2mat(header{:,6})) max(cell2mat(header{:,6}))];

if flag == 1, warning('Missing origin times.'); end
if flag == 2, warning('Missing S-wave times.'); end

%% VELOCITY MODELS: ORIGINAL, INVERSION, and PROPAGATION
% Save x,y,z in degrees switching as longitude comes second
% Find distance and azimuth to change in meters

wgs84                   =   wgs84Ellipsoid("m");
[d,az]                  =...
    distance(origin(1),origin(2),ending(1),ending(2),wgs84);

% convert az from degrees to radians for sin/cos
azr                     =   deg2rad(az);
dist_x                  =   d.*sin(azr);
dist_y                  =   d.*cos(azr);

% Coordinates of inversion points in meters
xM                      =   linspace(0,dist_x,nLong).';
yM                      =   linspace(0,dist_y,nLat).';
zM                      =   linspace(origin(3),ending(3),nzc).';
modvXYZ                 =   Murat_unfoldXYZ(xM,yM,zM);

Murat.input.x           =   linspace(origin(2),ending(2),nLong).';
Murat.input.y           =   linspace(origin(1),ending(1),nLat).';
Murat.input.z           =   linspace(origin(3),ending(3),nzc).';
Murat.input.gridStepX   =   xM(2) - xM(1);
Murat.input.gridStepY   =   yM(2) - yM(1);

if availableVelocity ==  0
    
    qLat                =   mean([origin(1),ending(1)]);
    qLon                =   mean([origin(2),ending(2)]);


    gridPropagation.x   =   xM';
    gridPropagation.y   =   yM';
    gridPropagation.z   =   zM';
    
    if isequal(Murat.input.namev,"LITHO1.0.nc")
        modvOriginal        =   Murat_readLithos1(velocityModel,qLat,qLon);
    else
        modvOriginal        =   load(velocityModel);
    end

    [modv,pvel,modvPlot]=...
        Murat_modv1D(modvXYZ,modvOriginal,PorS,origin);
    
    Murat.input.modv    =   modv;
    Murat.input.modvp   =   modv;
    Murat.input.modvPlot=   modvPlot;
    
    
elseif availableVelocity ==  1
    modvOriginal        =   load(velocityModel);

    [modvP,pvel,modvEqS,modvPlo]=   Murat_modv3D(FPath,FLabel,...
        modvOriginal,origin,mLat,mLon,nLat,nLong,nzc);
    
    gridPropagation.x   =   unique(modvP(:,1));
    gridPropagation.y   =   unique(modvP(:,2));
    gridPropagation.z   =   sort(unique(modvP(:,3)),'descend');
    
    Murat.input.modv    =   modvEqS;
    Murat.input.modvp   =   modvP;    
    Murat.input.modvPlot=   modvPlo;
end

Murat.input.gridPropagation     =   gridPropagation;
Murat.input.pvel                =   pvel;
Murat.input.header              =   header;
Murat.input.DDcoordinates       =   Murat_DDcoordinates(origin,ending,...
    nLat,nLong,nzc);
Murat.input.sacHeader           =   sacHeader;
outDirTXT                       =...
    fullfile(FLabel,'Tests','DataHeaders.xlsx');
writetable(header,outDirTXT);
end

%% Local helper (optimized)
function [listWithFolder, listNoFolder] = createsList(directory)
d                       = dir(directory);
if isempty(d)
    listWithFolder      = {};
    listNoFolder        = {};
    return
end
% remove hidden files (names starting with '.')
names                   = {d.name}.';
folders                 = {d.folder}.';
mask                    = ~startsWith(names, '.');
names                   = names(mask);
folders                 = folders(mask);
listWithFolder          = fullfile(folders, names);
listNoFolder            = names;
end