% ADDITIONAL input variables that are not set by the user.
function Murat                  =   Murat_checks(Murat,mag,timespan)

% INPUTS
Label                           =   ['./' Murat.input.label];
dataDirectory                   =   ['./' Murat.input.dataDirectory];

PTime                           =   ['SAChdr.times.' Murat.input.PTime];


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

origin                          =   Murat.input.origin;
ending                          =   Murat.input.end;
nLat                            =   Murat.input.gridLat;
nLong                           =   Murat.input.gridLong;
nzc                             =   Murat.input.gridZ;
velocityModel                   =   ['velocity_models/',Murat.input.namev];
sections                        =   Murat.input.sections;

if exist('./temp','dir')==7    
    delete('./temp/*')
else
    mkdir('./temp')
end

% Checking data
[Murat.input.listSac,~]         =	createsList(Murat,mag,timespan);
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
Murat.input.x                   =   linspace(origin(2),ending(2),nLong);
Murat.input.y                   =   linspace(origin(1),ending(1),nLat);
Murat.input.z                   =   linspace(origin(3),ending(3),nzc);

% Find distance and azimuth to change in meters - requires longitude first
dist_xdeg                       =   ending(2)-origin(2);
dist_ydeg                       =   ending(1)-origin(1);
dist_x                          =   deg2km(dist_xdeg)*1000;
dist_y                          =   deg2km(dist_ydeg)*1000;

% Coordinates of inversion points in meters
xM                         =   linspace(0,dist_x,nLong);
yM                         =   linspace(0,dist_y,nLat);
zM                         =   linspace(origin(3),ending(3),nzc);

Murat.input.gridStepX           =   xM(2)-xM(1);
Murat.input.gridStepY           =   yM(2)-yM(1);

if Murat.input.availableVelocity ==  0
    % The velocity model is a false 3D, it is used as
    % both inversion and propagation model
    modv                        =...
        Murat_unfoldXYZ(xM',yM',zM');
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
    gridPropagation.x           =   xM';
    gridPropagation.y           =   yM';
    gridPropagation.z           =   zM';
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
    % Original 3D velocity model from text file, in Lat/Log/Depth/V.
    modv_m                      =   load(velocityModel);
    
    modv_m1                     =   modv_m;

    modvDEG                     =...
        [modv_m1(:,1)-origin(1) modv_m1(:,2)-origin(2)...
        modv_m1(:,3:4)];

    % Switch as [lat,long] is [y,x] - the velocity in Murat format
    modv_o                      =   sortrows(...
        [deg2km(modvDEG(:,2))*1000 deg2km(modvDEG(:,1))*1000 ...
        modvDEG(:,3:4)],[1,2,-3]);

    % Nodes of the original velocity model
    xD                          =   unique(modv_o(:,1));
    yD                          =   unique(modv_o(:,2));
    zD                          =   sort(unique(modv_o(:,3)),'descend');
    
    % Uses meshgrid - fold switches again for interpolation and plot
    [XD,YD,ZD,V]                =   Murat_fold(xD,yD,zD,modv_o(:,4));
    
    % Interpolated axes for inversion model
    [X,Y,Z]                     =   meshgrid(xM,yM,zM);
    
    % Interpolated 3D inversion model
    mV                          =   griddata(XD,YD,ZD,V,X,Y,Z);
    
    %In case limits outside of the grid: interpolate
    if find(isnan(mV))
        mV                      =   inpaintn(mV);
    end
    Murat.input.modvPlot        =   mV;
    
    %Create the classical velocity model for inversion - again Murat format
    Murat.input.modv            =   Murat_unfold(X,Y,Z,mV);
    
    %Interpolated model for ray-tracing - half grid step for original
    resol2x                     =   abs(xM(2)-xM(1))/2;
    resol2y                     =   abs(yM(2)-yM(1))/2;
    resol2z                     =   abs(zM(2)-zM(1))/2;
    
    % Interpolated vectors
    xp                          =   xM(1):resol2x:xM(end);
    yp                          =   yM(1):resol2y:yM(end);
    zp                          =   zM(1):-resol2z:zM(end);
     
    % Meshes for interpolation
    [Xq,Yq,Zq]                  =   meshgrid(xp',yp',zp');
    
    % pvel is the propagation velocity model
    pvel                        =   griddata(XD,YD,ZD,V,Xq,Yq,Zq);
    if find(isnan(pvel))
        pvel                    =   inpaintn(pvel);
    end
    
    % Standard unpacking of the velocity model for propagation
    modvp                       =   Murat_unfold(Xq,Yq,Zq,pvel);
    Murat.input.modvp           =   modvp(:,1:4);
    
    % Nodes of the propagation model
    gridPropagation.x           =   unique(modvp(:,1));
    gridPropagation.y           =   unique(modvp(:,2));
    gridPropagation.z           =   sort(unique(modvp(:,3)),'descend');
    Murat.input.gridPropagation =   gridPropagation;
   
    %% plot models to check if they are correct
    
%     % check original model
%     Murat_test_vel_model('original velocity model',modv_m,sections)
%     % check centered model
%     Murat_test_vel_model('centered velocity model',modv_o,sections,origin)
%     %check interpolated version
%     Murat_test_vel_model('interpolated velocity model',mV,sections,origin,X,Y,Z)
%     %check interpolated version
%     Murat_test_vel_model('propagation velocity model',pvel,sections,origin,Xq,Yq,Zq)
    
end

Murat.input.pvel                =   pvel;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [listWithFolder,listNoFolder]...
    =   createsList(Murat,mag,timespan)
% CREATES a list of visible files in a folder, outputs both with and
% without folder
%
if isempty(mag) && isempty(timespan)
    list                            =   dir(Murat.input.dataDirectory);
    list                            =   list(~startsWith({list.name}, '.'));
    
    listWithFolder                  =	fullfile({list.folder},{list.name})';
    listNoFolder                    =   {list.name}';
    
else
    %% Miriam's option to only read events of certain magnitude
    
    DFolder                     = Murat.input.dataDirectory;
    load('chosen_events.mat')
    orig_data = length(TZ_events.year);
    
    % delete traces outside main area of interest
    to_del = zeros(length(TZ_events.year),1);
    for i=1:length(TZ_events.year)
        if TZ_events.lat(i) > -2.45 || TZ_events.lat(i) < -2.95
            to_del(i) = 1;
        end
        if TZ_events.lon(i) < 35.75  || TZ_events.lon(i) > 36.25
            to_del(i) = 1;
        end
    end
    del_row = find(to_del==1);
    TZ_events(del_row,:)=[];
    disp(['removed ' , num2str(length(find(to_del==1))), ' events because they are far from network'])
    % delete traces below magnitude threshold
    find_mags = find(TZ_events.mag<mag);
    TZ_events(find_mags,:) = [];
    disp(['removed ' , num2str(length(find_mags)), ' events because they are below the magnitude threshold'])
    % keep traces inside chosen time frame
    ex_dates=cellstr(TZ_events.of);
    dates = cellfun(@(x) x(1:6),ex_dates,'un',0);
    find_dates = contains(dates,timespan);
    TZ_events(~find_dates,:) = [];
    disp(['selected ' , num2str(length(find(find_dates==1))), ' events inside the set timeframe'])
    % delete traces below certain depth
    find_depth = find(TZ_events.depth>25);
    TZ_events(find_depth,:) = [];
    disp(['removed ' , num2str(length(find_depth)), ' events because they are below 25 km depth'])
    % delete events which have only few traces
    find_ns = find(TZ_events.ns<10);
    TZ_events(find_ns,:) = [];
    disp(['removed ' , num2str(length(find_ns)), ' events because they were recorded by less than 10 components'])
    
    disp(['keeping ',num2str(length(TZ_events.lat)), 'of ', num2str(orig_data),' events'])
    
    new_filenames = [];
    new_filedir = [];
    
    for ii = 1:length(TZ_events.year)
        tmp_filename = sprintf('%4d%02d%02d%02d%02d%02d',TZ_events.year(ii),TZ_events.month(ii),...
            TZ_events.day(ii),TZ_events.hour(ii),TZ_events.min(ii),floor(TZ_events.secs(ii)));
        tmp_lst = dir([DFolder,tmp_filename,'*']);
        tmp_filenames = {tmp_lst.name}';
        tmp_filedir = {tmp_lst.folder}';
        new_filenames = [new_filenames;tmp_filenames];
        new_filedir = [new_filedir;tmp_filedir];
    end
    
    listWithFolder = fullfile(new_filedir,new_filenames);
    listNoFolder = new_filenames;
    
end
