% PROGRAM FOR: 3D direct- and 2D peak delay and coda-wave (with and without
% kernels) attenuation tomography
%
% Author:  L. De Siena, January 2018 
% 
% Installation
% ------------
% 
% SYSTEM: The program works on a Macbook pro with High Sierra, Matlab R2017a.
% Necessary Toolboxes: Signal Processing and Mapping (for latitute-longitude coordinates - Romania example).
% 
% Two sample datasets (Mount St. Helens and Romania) can be downloaded at https://doi.pangaea.de/10.1594/PANGAEA.893893 to test the code. A third sample input (Pollino) is provided and the related dataset may be requested to Luca De Siena.
% 
% INSTRUCTIONS:
% The current version works following these steps:
% 
% 1. Download the package at https://github.com/LucaDeSiena/MuRAT.
% 
% 2. Download the two sample datasets at
% https://doi.pangaea.de/10.1594/PANGAEA.893893.
% Unzip the MSH (Mount St. Helens) and Romania datasets and put
% the folders in the Murat-master folder.
% 
% 3. Build your own input file (.xsl) - each field is described in the
% attached README file and in the INPUT sections of thi code.
% When building your example, use one of the Input
% files as template (MSH for 3D and Romania for 2D). Always start with a
% 2D analysis (“pa=1” or “pa=2”). “pa=3” only works with a suitable
% velocity model.
% 
% 4. Run MuRAT2.m.
% 
% 5. When asked insert the name of the input file desired
% (either Input_MSH.m or Input_Romania.m).
% 
% 6. After the L curves are produced, write the smoothing parameter.
%  
%
% Theory:
% ----
% Peak-delays and 2D coda attenuation: Takahashi et al. 2007 (JGR);
% Calvet et al. 2013b (GJI); De Siena et al. 2016 (EPSL);
% Borleanu et al. 2017 (Tectonophysics)
%
% Kernel-based 2D coda attenuation: Del Pezzo et al. 2016 (GJI);
% De Siena et al. 2017 (GRL)
%
% Direct-wave coda-normalized attenuation: Del Pezzo et al. 2006 (PEPI);
% De Siena et al. 2010 (JGR); De Siena et al. 2014a (MuRAT, JVGR) and
% De Siena et al. 2014b (JGR)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUTS AND CHECKS

unzip('Utilities_Matlab.zip','Utilities_Matlab')
addpath('./Utilities_Matlab')
clear
close all
clc

% Import data from spreadsheet
%
%    Sample Workbook: ./Input_file_MSH.xlsx or Input_file_Romania.xlsx

[file,path,indx] = uigetfile('*.*');
[~, ~, raw] = xlsread(fullfile(path, file));
raw = raw(end,:);
stringVectors = string(raw(:,[2,3,4,5,6,7,19,20,21,25,26,27,35]));
stringVectors(ismissing(stringVectors)) = '';
raw = raw(:,[1,8,9,10,11,12,13,14,15,16,17,18,22,23,24,28,29,30,31,32,...
    33,34,36,37,38,39,40,41,42,43,44,45]);

% Create output variable
data = reshape([raw{:}],size(raw));

% Create table
InputfileMSH = table;

% Which analysis do you want to perform?
% Pick delay and Qc without kernels: pa=1
% Pick delay and Qc with kernels: pa=2
% Pick delay, kernel-Qc and P/S wave attenuation with the CN method: pa=3
pa = data(:,1);
Murat.analysis = pa;

%% INPUT PATHS AND FORMATS
% Creates folders and paths to store results and figures
Murat.paths.label = stringVectors{:,1};
Murat.paths.workingdir = stringVectors{:,2};

% Creating a list of files in folder
% Folder containing data
DFolder = stringVectors{:,3};
Murat.paths.datadir =DFolder;

% Output figure format 
Murat.figures.format = stringVectors{:,4};

% Figure visibility - pick 'on' or 'off'
Murat.figures.visibility = stringVectors{:,5};

%set here the planes you want to inspect - only for pa=3
if pa==3
    sections = str2num(stringVectors{:,6});
    Murat.figures.sections=sections;
end

% How big are markers for stations and events
Murat.figures.sizeMarker=60;

% Check that user has the Mapping Toolbox installed.
hasMT = license('test', 'map_toolbox');
Murat.figures.hasMT=hasMT;
if ~hasMT
  % User does not have the toolbox installed.
  message =...
      sprintf('No Mapping Toolbox - Maps will not be geo-localised');
end

%% INPUT DATA
% Choose between P- (2) and S- (3) waves for peak delay
% RULE: First column in time files always contains origin of events
% RULE: Second column in time files always contains P-wave phase
% RULE: Third column in time files always contains S-wave phase
Murat.data.PorS = data(:,2);

% Work with 1 vertical (1) or 2 horizontal (2) recordings, or with the
% three components(3)
Murat.data.components =  data(:,3);

% To see spectrogram (set>0) of the seismogram: set seisi='number of
% seismogram to see spectrogram'
Murat.data.spectrogram = data(:,4);

% Central frequency - set it according to your spectrograms
cf = data(:,5);
Murat.data.centralFrequency = cf;

% Maximum window to pick pick-delays in seconds
Murat.data.maximumPD = data(:,6);

% Minimum peak delay considering scattering in the area and frequency
Murat.data.minimumPD = data(:,7);

% Lapse time for the start of the window used to measure and calculate the
% normalization energy
Murat.data.startLT = data(:,8);

% Total coda window length for Qc and normalization
Murat.data.codaWindow = data(:,9);

% Parameter for smoothing - must be > 2
Murat.data.smoothing = data(:,10);

% The sped coefficient sets the spectral energy decay of the coda
% wavefield
Murat.data.spectralDecay = data(:,11);

% Number of nodes along WE and SN
% If the origin time is unknown, you can set a theoretichal velocity for
% the whole area and evaluate it from picking. It must be the velocity of
% the phase you are mapping. Velocity is in km/s.
Murat.data.averageVelocity = data(:,12);

% Name of the SAC variables where zero time and P/S pickings are saved
Murat.data.originTime = stringVectors{:,7};
Murat.data.PTime = stringVectors{:,8};
Murat.data.STime = stringVectors{:,9};

%start time to measure noise - before P-arrival in seconds - for pa=3
if pa==3
    % Length of the window used to measure noise and direct wave intensity
    Murat.data.bodyWindow = data(:,13);

    % Start of the window used to measure noise
    Murat.data.startNoise = data(:,14);
    
    % The minimum coda-to-noise energy ratio for the weighted inversion
    Murat.data.tresholdNoise = 10;

end

%% INPUT GEOMETRY
% Import event origin time and coords of event and station from files
% (Import=1).
% Ideal is to have both in the SAC header (Import=2) and do it
% in lat/long.
evst = data(:,15);
Murat.geometry.import=evst;

% UTM coordinates of the origin of the model - this must be put at
% least one cell before the westernmost and southernmost
% earthquakes/station
origin = str2num(stringVectors{:,10});
Murat.geometry.origin=origin;

%Step of the grid and number of nodes for pd and Qc imaging
nxc = data(:,16);
Murat.geometry.gridX = nxc;

nyc = data(:,17);
Murat.geometry.gridY = nyc;

stepg = data(:,18);
Murat.geometry.gridStep=stepg;

degorutm = data(:,19);
Murat.geometry.degreesorutm=degorutm;

if evst==1
    
    % Opening the event file. Format is:
    % column (1) = twelve numbers for the origin time of the event
    % (date+time in seconds)
    % column (2) = UTM (WE) or latitude
    % column (3) = UTM (SN) or longitude
    % column (4) = Altitude above sea level in meters
    event = fopen(stringVectors{:,11});
    namee= textscan(event,'%s %f %f %f');
    Murat.geometry.nameeven=namee{1};
    Murat.geometry.even=[namee{2} namee{3} namee{4}];
    fclose(event);
    
    % Opening the station file. Format is:
    % column (1) = Name of station (3 characters)
    % column (2) = UTM (WE) or latitude
    % column (3) = UTM (SN) or longitude
    % column (4) = Altitude above sea level in meters
    station=fopen(stringVectors{:,12});
    names= textscan(event,'%s %f %f %f');
    Murat.geometry.namestation=names{1};
    Murat.geometry.station=[names{2} names{3} names{4}];
    fclose(station);
    
end

%pre-define 2D matrix in space
XY = zeros(floor(nxc)*floor(nyc),2);

% Some figures require moving the point of plotting to half the resolution
halfres=stepg/2;
Murat.geometry.x=origin(1)+halfres:stepg:stepg*(nxc-1)+origin(1)+halfres;
Murat.geometry.y=origin(2)+halfres:stepg:stepg*(nyc-1)+origin(2)+halfres;

index=0;
for i=1:nxc
    for j=1:nyc
        index=index+1;
        XY(index,1:2)=[origin(1)+stepg*i origin(2)+stepg*j];
    end
end
Murat.geometry.map=XY;
lxy=length(XY(:,1));

% Inputs necessary for the direct-wave CN attenuation tomography
if Murat.analysis==3
    
    %set createrays to 1 if you want to create files for rays
    %WARNING: This will create A BIG .mat file
    Murat.geometry.createrays = data(:,20);
    
    dtreshold = data(:,21);
    Murat.geometry.depthTreshold=dtreshold;
    
    %Rays measured in meters or degrees.    
    if degorutm==1
        Murat.geometry.unity = 1000;
    elseif degorutm==111
        Murat.geometry.unity = 1;
    end
    
    %1D (1) or 3D (3) velocity model
    dimV = data(:,22);

    % This works in [x,y,z], created for UTM WGS84 and origins to zero.
    % In the new version if a 3D velocity model is unavailable a false 3D
    % is created from iasp91
    if dimV==1
        
        modv1=load('iasp91.txt'); %1D velocity model iasp91
        dend=-34000; %select maximum depth range
        
        li=length(modv1(:,1));
        modv=zeros(lxy*li,5);
        index=0;
        for i=1:nxc
            for j=1:nyc
                index1=index;
                index=index+1;
                modv(index1*li+1:index1*li+li,1)=XY(index,1)*1000;
                modv(index1*li+1:index1*li+li,2)=XY(index,2)*1000;
                modv(index1*li+1:index1*li+li,3)=...
                    -modv1(1:li,1)*1000;
                modv(index1*li+1:index1*li+li,4)=...
                    modv1(1:li,PorS+1);
            end
        end
        
        modv(:,1)=(modv(:,1)-modv(1,1));
        modv(:,2)=(modv(:,2)-modv(1,2));
        modv(modv(:,3)<dend,:)=[];
        
    elseif dimV==3
    
        modv=load(stringVectors{:,13}); %3D velocity model from text file
        modv(:,5)=0;
        modv(:,1)=modv(:,1)+origin(1);
        modv(:,2)=modv(:,2)+origin(2);
        
    end
    
    Murat.geometry.modv=modv;
    
    chx=find(modv(:,1)~=modv(1,1),1);
    chy=find(modv(:,2)~=modv(1,2),1);
    resol2x = abs(modv(chx,1)-modv(1,1))/2;
    resol2y = abs(modv(chy,2)-modv(1,2))/2;
    resol2z = abs(modv(2,3)-modv(1,3))/2;
    
    
    %Steps of the velocity models
    passox=max(modv(:,1))-min(modv(:,1));
    passoy=max(modv(:,2))-min(modv(:,2));
    passoz=max(modv(:,3))-min(modv(:,3));
    
    passo=[passox passoy passoz];
    resol=[resol2x resol2y resol2z];
    resol2=min(resol);
    Murat.geometry.resolutionMin=resol2;
    
    %Creates grid for direct waves and check for zeroes
    %Regular step of the gridD for interpolation - half of step of modv
    
    ixD=floor(passox/resol2x);%numer of x,given the step of gridD
    iyD=floor(passoy/resol2y);%numer of y,given the step of gridD
    izD=floor(passoz/resol2z);%numer of depths,given the step of gridD
    
    gridD=zeros(3,max(passo./resol));
    gridD(1,1:ixD)=origin(1):resol2x:origin(1)+passox-resol2x;
    gridD(2,1:iyD)=origin(2):resol2y:origin(2)+passoy-resol2y;
    gridD(3,1:izD)=-origin(3):resol2z:-origin(3)+passoz-resol2z;
    Murat.geometry.gridD=gridD;
    % gridD dimensions
    pvel = zeros(ixD,iyD,izD);
    
    %NUMBER OF X, Y, AND Z LAYERS
    for k=1:izD
        index=0;
        for i=1:ixD
            for j=1:iyD
                index=index+1;
                pvel(i,j,k) = modv(index,4);
            end
        end
    end
    Murat.geometry.pvel=pvel;
end

% Name of the file containing the seismograms and rays
% creates a list of the waveforms and calculates their number.
% RULE: This is necessary for the option Murat.geometry.import = 1, where
% event and station info are in two separate files. In this case,
% the SAC file-name must have a specific format:
% the first 12 characters are the origin time of the event, as per
% event-location file). Characters 18:20 are the name of station as per
% station file.
%
% For Murat.geometry.import = 2, there is no requirement on SAC filenames.
% In the case of pa=3 the names of rays created will be identical to files.

% Get name of files in traces' folder
list = dir(DFolder);
filenames = {list.name}'; %create cell array of file names
filedir = {list.folder}'; %create cell array of file folder
listasac = filenames;
lls = length(listasac);
lung=zeros(lls,1);
lista = filenames; %This will be the name of the file listing rays

for i = 1:lls
    listasac{i,1} = cat(2,filedir{i},'/',filenames{i});
    %Here it takes the event/station name. It is necessary to adapt the
    %numbers to where even name and station name are
    if evst==1
        li = lista{i,1};
        li1 = cat(2,li(1:12),li(18:20));
        lista{i,1} = li1;
    end
end

Murat.paths.listasac=listasac;
Murat.paths.lista=lista;

%% INPUT INVERSION
% Size anomaly for testing: twice(2). Might be set to four(4). The input of
% the checkerboard must be always checked visually at the end of
% the process    
Murat.inversion.sizeCheck = data(:,23);

% Values of attenuation for testing
Murat.inversion.highCheck = data(:,24);
Murat.inversion.lowCheck = data(:,25);

% Seconds of time-windows for non-linear inversion and corresponding
% number, set nonlinear=1 to activate;
Murat.inversion.nonlinear = data(:,26);

% Uncertainty on Qc estimation
if Murat.inversion.nonlinear==0
    % Minimum R-squared for Qc fitting
    Murat.inversion.fitT = data(:,27);
elseif Murat.inversion.nonlinear==1
    % Length of smaller time windows to compute compute coda intensity
    Murat.inversion.fitL = data(:,28);
    Murat.inversion.fitT = data(:,29);
    % Number of time windows to compute coda intensity
    Murat.inversion.fitN=Murat.data.codaWindow/Murat.inversion.fitL;
    %Grid search - set the minimum, maximum and spacing to search in the
    %parameter space
    m1min = data(:,30);
    Murat.inversion.minimum=m1min; % minimum inverse Qc allowed
    m1max = data(:,31);
    Murat.inversion.maximum=m1max;% minimum inverse Qc allowed
    total = data(:,32);
    Murat.inversion.total = total;% total number of Qc in the interval
    Murat.inversion.fit = (m1min + (m1max-m1min)/(total-1)*(0:total-1))';
end

if exist(Murat.paths.label,'dir')~=7
    mkdir(Murat.paths.label)
end

clearvars -except Murat
save('Murat.mat','Murat');

% The following prompts will measure peak-delays (pd) ad Qc depending on
% the input files.
%% RAY TRACING METHOD: RAY BENDING - for evst=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% With option evst=1, we use two external files for event and station
% locations. The SAC file name must have a specified format for recognition
disp('Geometry Section')

load Murat.mat

pa=Murat.analysis;

%PATHS and FIGURES
FPath=Murat.paths.workingdir;
FLabel=Murat.paths.label;
visib=Murat.figures.visibility;
sz=Murat.figures.sizeMarker;
fformat=Murat.figures.format;

%DATA
evst=Murat.geometry.import;
lista=Murat.paths.lista;
lls=length(lista);

%GEOMETRY
XY=Murat.geometry.map;
lxy=length(XY(:,1));
origin=Murat.geometry.origin;
nxc=Murat.geometry.gridX;
nyc=Murat.geometry.gridY;
x=Murat.geometry.x;
y=Murat.geometry.y;

if pa==3
    gridD=Murat.geometry.gridD;
    pvel=Murat.geometry.pvel;
    modv=Murat.geometry.modv;
    um=Murat.geometry.unity;
    if evst ==1
        numev=length(Murat.geometry.nameeven);
        numst=length(Murat.geometry.namestation);
        createrays=Murat.geometry.createrays;
        nameeven=Murat.geometry.nameeven;
        even=Murat.geometry.even;
        namesta=Murat.geometry.namestation;
        staz=Murat.geometry.station;
    end
end

if evst==1
    %     tic
    disp('Inversion Section for events and stations from files')
    
    %Predefine the inversion matrix for pd and Qc imaging
    Ac=zeros(lls,lxy);
    Apd=zeros(lls,lxy);
    D=zeros(lls,1);
    evestaz=zeros(lls,4);
    
    % This loop creates both rays and elements of Apd and Ac - for peak
    % delays and Qc imaging in case of files of events and stations
    % outside SAC haeder
    
    indexray=0; % Increases every time a ray is created.
    
    %Sets figure for rays - it will show rays in the reference system of
    %the pre-defined grid. If the Mapping Toolbox is available and the
    %coordinates are in latitude and longitude, the figure will show coast
    %lines
    rays=figure('Name','Rays','NumberTitle','off','visible',visib);
    
    if pa==3
        lunparz = zeros(100,lls);
        blocchi = zeros(100,lls);
        sb = zeros(100,lls);
        luntot = zeros(lls,1);
        
        if createrays ==1
            ma1=zeros(100,5,lls);
        end
    end
    % Loop over events in even.txt file
    for ii=1:numev
        ray(1,1)=even(ii,1);
        ray(2,1)=even(ii,2);
        ray(3,1)=-even(ii,3);
        namee1=nameeven{ii};
        namee=namee1(1:12);
        
        % Loop over stations in staz.txt file
        for ir=1:numst
            namest1=namesta{ir};
            namest=namest1(1:3);
            c1=cat(2,namee,namest);
            
            % Pattern recognition from the name of the SAC file
            if find(strncmp(lista(:,1),c1,length(c1)))>0
                indexray=indexray+1;
                if isequal(mod(indexray,100),0)
                    disp(['Ray number is ', num2str(indexray)])
                end
                ray(1,2) = staz(ir,1);
                ray(2,2) = staz(ir,2);
                ray(3,2) = -staz(ir,3);
                xx=XY(:,1);
                yy=XY(:,2);
                miix=min(ray(1,1),ray(1,2));
                maax=max(ray(1,1),ray(1,2));
                miiy=min(ray(2,1),ray(2,2));
                maay=max(ray(2,1),ray(2,2));
                fxy=find(xx>miix & xx<maax & yy>miiy & yy<maay);
                
                Apd(indexray,fxy)=1;% Creates pd matrix
                sst = [ray(1,1) ray(2,1) ray(1,2) ray(2,2)];
                evestaz(indexray,1:4)=sst;% Creates file of selected coords
                
                % Epicantral distance
                D(indexray,1)=...
                    sqrt((sst(3)-sst(1))^2+(sst(4)-sst(2))^2)/1000;
                
                if pa==1
                    
                    % Same pd and Qc inversion matrix
                    Ac=Apd;
                    
                elseif pa>1
                    
                    % Operations to calculate normalised kernels
                    D1=sqrt((sst(3)-sst(1))^2+(sst(4)-sst(2))^2);
                    deltaxy=0.2;
                    F1=1/2/pi/deltaxy^2/D1^2;
                    F2=1/deltaxy^2/D1^2;
                    F3=F1*(0.5*exp(-abs((xx-(sst(1)+sst(3))/2).^2*F2/2+...
                        (yy-(sst(2)+sst(4))/2).^2*F2/0.5))+...
                        exp(-abs((xx-sst(1)).^2*F2/2+...
                        (yy-sst(2)).^2*F2/2))+...
                        exp(-abs((xx-sst(3)).^2*F2/2+...
                        (yy-sst(4)).^2*F2/2)));
                    if find(F3)>0
                        F=F3/sum(F3);
                    else
                        F=F3;
                    end
                    no=F<0.0001;
                    F(no)=0;
                    
                    %Inversion matrix for Qc
                    Ac(indexray,:)=F;
                    
                end
                
                % In the case a velocity model is available
                if pa==3
                    
                    % Function for ray-bending
                    rma=tracing(ray,gridD,pvel);
                    
                    % Set this to create reay-files
                    if createrays ==1
                        lrma=length(rma(:,1));
                        ma1(1:lrma,:,indexray) = rma;
                    end
                    
                    %Creates figure with rays
                    hold on
                    subplot(2,2,1)
                    plot(rma(:,2),rma(:,3),'k');
                    xlabel('WE','FontSize',12,'FontWeight','bold',...
                        'Color','k')
                    ylabel('SN','FontSize',12,'FontWeight','bold',...
                    'Color','k')
                    hold off
                    hold on
                    subplot(2,2,2)
                    plot(rma(:,4),rma(:,3),'k');
                    xlabel('Depth UTM (WGS84)','FontSize',12,...
                        'FontWeight','bold','Color','k')
                    ylabel('SN','FontSize',12,'FontWeight','bold',...
                        'Color','k')
                    hold off
                    hold on
                    subplot(2,2,3)
                    plot(rma(:,2),-rma(:,4),'k');
                    xlabel('WE','FontSize',12,'FontWeight','bold',...
                        'Color','k')
                    ylabel('Depth UTM (WGS84)','FontSize',12,...
                        'FontWeight','bold','Color','k')
                    hold off
                    
                    %======================================================
                    % We calculate segments in a grid, defined from the
                    % velocity model. The rays are in the reference system
                    % of the velocity model, first point at the source,
                    % and with positive depth. The velocity model is a
                    % x,y,z grid, with the same steps in the three
                    % directions. The outputs are necessary to invert for
                    % direct-wave attenuation - Analysis = 3
                    %======================================================
                    
                    [lunpar, blocch, lunto, s, modv] =...
                        segments_single(modv,origin(3),um,rma);
                    lunparz(1:length(lunpar),indexray)=lunpar;
                    blocchi(1:length(blocch),indexray)=blocch;
                    luntot(indexray,1)=lunto;
                    sb(1:length(s),indexray)=s;
                    
                end
            end
        end
    end
    Murat.inversion.APeakDelay=Apd;
    Murat.inversion.AQCoda=Ac;
    Murat.inversion.epicentralDistance=D;
    
    if pa<3
        
        %Create figure with 2D rays
        for nn=1:lls
            hold on
            plot([evestaz(nn,1) evestaz(nn,3)],...
                [evestaz(nn,2) evestaz(nn,4)],'k-')
        end
        hold on
        scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
        hold on
        scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
        hold off
        
        grid on
        ax = gca;
        ax.GridLineStyle = '-';
        ax.GridColor = 'k';
        ax.GridAlpha = 1;
        ax.LineWidth = 1;
        
    elseif pa==3
        
        hold on
        
        %Create figure with 3D rays
        subplot(2,2,1)
        scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold on
        scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold off
        hold on
        subplot(2,2,2)
        scatter(-even(:,3),even(:,2),sz,'c','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold on
        scatter(-staz(:,3),staz(:,2),sz,'^','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold off
        hold on
        subplot(2,2,3)
        scatter(even(:,1),even(:,3),sz,'c','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold on
        scatter(staz(:,1),staz(:,3),sz,'^','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold off
        
        % Only solve for the blocks crossed by at least 1 ray
        fover = modv(:,5)>=1;
        inblocchi = modv(fover,:);
        inblocchi(:,5)=find(fover);
        
        %INVERSION MATRIX for direct waves: A
        A = zeros(length(Ac(:,1)),length(inblocchi(:,1)));
        for j = 1:length(sb(1,:))
            for i = 2:length(sb(:,1))
                bl = blocchi(i,j);
                slo = sb(i,j);
                l = lunparz(i,j);
                fbl = find(inblocchi(:,5)==bl);
                if isempty(fbl)==0
                    % check the inversion matrix element IN SECONDS
                    A(j,fbl)=-l*slo;
                end
            end
        end
        
        Murat.inversion.AQ=A;
        Murat.inversion.partialLengths=lunparz;
        Murat.inversion.totalLengths=luntot;
        Murat.inversion.hitBlocks=inblocchi;
        Murat.inversion.slownessBlocks=sb;
        
        if createrays == 1
            Murat.results.rays=ma1;
        end
    end
    
    % Saves the figure of RAYS
    FName = 'Rays';
    saveas(rays,fullfile(FPath, FLabel, FName), fformat);
    
    if pa>1
        
        % Kernel sensitivity figure
        Qcsen=figure('Name','Qc sensitivity, first source-station pair',...
            'NumberTitle','off','visible',visib);
        Qcss=Ac(1,:);
        Qcs=zeros(nxc,nyc);
        
        index=0;
        for i=1:length(x)
            for j=1:length(y)
                index=index+1;
                Qcs(i,j)=Qcss(index);
            end
        end
        
        [X,Y]=meshgrid(x,y);
        contourf(X,Y,Qcs')
        
        c = colorbar;
        c.Label.String = 'Qc';
        grid on
        ax = gca;
        ax.GridLineStyle = '-';
        ax.GridColor = 'k';
        ax.GridAlpha = 1;
        ax.LineWidth = 1;
        
        FName = 'Qc_sensitivity';
        saveas(Qcsen,fullfile(FPath, FLabel, FName), fformat);
    end
end

clearvars -except Murat
save('Murat.mat','Murat');
% rtoc=toc;
%% Seismic attributes for peak delay, Qc and Qp,s imaging
clear
% tic
disp('Data Section')

load Murat.mat

pa=Murat.analysis;

%PATHS and FIGURES
FPath=Murat.paths.workingdir;
FLabel=Murat.paths.label;
visib=Murat.figures.visibility;
sz=Murat.figures.sizeMarker;
fformat=Murat.figures.format;
hasMT=Murat.figures.hasMT;

%DATA
evst=Murat.geometry.import;
listasac=Murat.paths.listasac;
lls=length(listasac);
tCm=Murat.data.startLT;
tWm=Murat.data.codaWindow;
f0=Murat.data.originTime;
fP=Murat.data.PTime;
fS=Murat.data.STime;
seisi=Murat.data.spectrogram;
cf=Murat.data.centralFrequency;
mintpde=Murat.data.minimumPD;
maxtpde=Murat.data.maximumPD;
nf=Murat.data.smoothing;
sped=Murat.data.spectralDecay;
vth=Murat.data.averageVelocity;


%GEOMETRY
PorS=Murat.data.PorS;
compon=Murat.data.components;
XY=Murat.geometry.map;
lxy=length(XY(:,1));
origin=Murat.geometry.origin;
nxc=Murat.geometry.gridX;
nyc=Murat.geometry.gridY;
stepg=Murat.geometry.gridStep;
x=Murat.geometry.x;
y=Murat.geometry.y;
degorutm=Murat.geometry.degreesorutm;


% INVERSION
nW=Murat.inversion.fitL;
ntW=Murat.inversion.fitN;
RZZ2=Murat.inversion.fitT;
m1a=Murat.inversion.fit;
nonlinear=Murat.inversion.nonlinear;
L1=Murat.inversion.total;
        
if pa==3
    fin1=Murat.data.bodyWindow;
    sn=Murat.data.startNoise;
    tresholdnoise=Murat.data.tresholdNoise;
    gridD=Murat.geometry.gridD;
    pvel=Murat.geometry.pvel;
    modv=Murat.geometry.modv;
    um=Murat.geometry.unity;
        
    if evst ==1
        createrays=Murat.geometry.createrays;
        even=Murat.geometry.even;
        staz=Murat.geometry.station;
        lunparz = Murat.inversion.partialLengths;
        luntot  = Murat.inversion.totalLengths;
        sb = Murat.inversion.slownessBlocks;
        A = Murat.inversion.AQ;
    end
end

%==========================================================================
% A loop to measure the Qc, peak-delay, and energy ratios for
% each seismic trace located in
% the folder "traces" with the coda-normalization method.
% The seismograms may be the output of both single and three-component
% seismic stations. In case of more than one component the rays must be
% ordered in the folder starting with the WE component, then SN and finally
% the Z component for each ray.
% 
% The program accepts SAC files and uses the information in the header
% The header must include peaking of the direct phase of interest, marker
% "a" for P and "t0" for S. If evst=2 you can get event and station info
% directly from the header. This will create the rays and related files
% automatically

if pa==3
    signal = zeros(lls,1); % The direct wave energies
    coda = zeros(lls,1); % The coda wave energies
    rapsp = zeros(lls,1); % The spectral ratios
    rapspcn = zeros(lls,1); % The coda versus noise ratios
end

tempi= zeros(lls,4); %Variables containing time origin (column 1), 
Qm = zeros(lls,1); %Qc with linearised or non-linear theory
RZZ = zeros(lls,1); %Correlation coefficient with respect to linear
ttheory=zeros(lls,1);
peakd=zeros(lls,1);
constQmean=zeros(2,2);
tlapse=(tCm+nW/2:nW:tCm+tWm-nW/2)';

index=0;
indexlonger=0;
    
if evst==1
    Ac=Murat.inversion.AQCoda;
    Apd=Murat.inversion.APeakDelay;
    D=Murat.inversion.epicentralDistance;
    
    if pa==3
        lunparz = Murat.inversion.partialLengths;
        blocchi = Murat.inversion.hitBlocks;
        sb = Murat.inversion.slownessBlocks;
        luntot = Murat.inversion.totalLengths;
    end
        
elseif evst==2
    disp('& Inversion Section for evst = 2')
    %Predefine the inversion matrix for peak-delay and coda-Q imaging
    even=zeros(lls,3);
    staz=zeros(lls,3);
    Ac=zeros(lls,lxy);
    Apd=zeros(lls,lxy);
    D=zeros(lls,1);
    evestaz=zeros(lls,1);
    rays=figure('Name','Rays','NumberTitle','off','visible',visib);
        
    if pa==3
        lunparz = zeros(100,lls);
        blocchi = zeros(100,lls);
        sb = zeros(100,lls);
        luntot = zeros(lls,1);
        if createrays ==1
            ma1=zeros(100,5,lls);
        end
    end
end

for i=1:lls
    index = index+1;
    if isequal(mod(index,100),0)
        disp(['Waveform number is ', num2str(index)])
    end
    
    %Read seismogram and get peaking/event/station info
    [tempis,sisma,SAChdr] = fget_sac(listasac{i}); % Converts SAC files
    
    tempi(i,2)=eval(fP);
    if tempi(i,2)==-12345
        continue
    end
    if isequal(fS,'[]')==0
        tempi(i,3)=eval(fS);
    end
    srate=1/SAChdr.times.delta; %sampling frequency
    
    if i==seisi %to see the seisi spectrogram
        spect=figure('Name','Spectrogram','NumberTitle','off',...
            'visible',visib);
        view(2)
        title('Spectrogram of the first seism','FontSize',12,...
            'FontWeight','bold','Color','k');
        spectrogram(sisma,50,30,50,srate,'yaxis')
        FName = 'Spectrogram';
        saveas(spect,fullfile(FPath, FLabel, FName), fformat);
    end
    
    % Filter
    Wn = ([cf-cf/3 cf+cf/3]/srate*2); %frequency band
    [z,p,k] = butter(4,Wn,'bandpass'); %butter filter
    % [z,p,k] = cheby1(4,4,Wn,'bandpass'); %chebyshev filter
    [sos,g] = zp2sos(z,p,k); % Convert to SOS form
    Hd = dfilt.df2tsos(sos,g); % Create a dfilt object   

    if evst==2
        even1=[SAChdr.event.evlo SAChdr.event.evla SAChdr.event.evdp];
        staz1=[SAChdr.station.stlo SAChdr.station.stla...
            SAChdr.station.stel];
        even(i,1:3)=even1; 
        staz(i,1:3)=staz1; 
        ray(1,1) = even1(1);
        ray(2,1) = even1(2);
        ray(3,1) = -even1(3);
        ray(1,2) = staz1(1);
        ray(2,2) = staz1(2);
        ray(3,2) = -staz1(3);
        xx=XY(:,1);
        yy=XY(:,2);    
        miix=min(even1(1),staz1(1));
        maax=max(even1(1),staz1(1));
        miiy=min(even1(2),staz1(2));
        maay=max(even1(2),staz1(2));
        fxy=find(xx>miix & xx<maax & yy>miiy & yy<maay);
        Apd(i,fxy)=1;
        sst = [even1(1) even1(2) staz1(1) staz1(2)];
        evestaz(i,1:4)=sst;
        D(i,1)=sqrt((sst(3)-sst(1))^2+(sst(4)-sst(2))^2);
        if pa==1
            
            % Same pd and Qc inversion matrix
            Ac=Apd;
            
        elseif pa>1
            
            % Operations to calculate normalised kernels
            D1=sqrt((sst(3)-sst(1))^2+(sst(4)-sst(2))^2);
            deltaxy=0.2;
            F1=1/2/pi/deltaxy^2/D1^2;
            F2=1/deltaxy^2/D1^2;
            F3=F1*(0.5*exp(-abs((xx-(sst(1)+sst(3))/2).^2*F2/2+...
                (yy-(sst(2)+sst(4))/2).^2*F2/0.5))+...
                exp(-abs((xx-sst(1)).^2*F2/2+...
                (yy-sst(2)).^2*F2/2))+...
                exp(-abs((xx-sst(3)).^2*F2/2+...
                (yy-sst(4)).^2*F2/2)));
            if find(F3)>0
                F=F3/sum(F3);
            else
                F=F3;
            end
            no=F<0.0001;
            F(no)=0;
            
            %Inversion matrix for Qc
            Ac(i,:)=F;
            
        end
        
        if pa==3
            rma=tracing(ray,gridD,pvel);
            if createrays ==1
                lrma=length(rma(:,1));
                ma1(1:lrma,:,i) = rma;
            end
                
            lrma=length(rma(:,1));
            
            %Creates figure with rays
            hold on
            subplot(2,2,1)
            plot(rma(:,2),rma(:,3),'k');
            xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
            ylabel('SN','FontSize',12,'FontWeight','bold','Color','k')
            hold off
            hold on
            subplot(2,2,2)
            plot(rma(:,4),rma(:,3),'k');
            xlabel('Depth UTM (WGS84)','FontSize',12,...
                'FontWeight','bold','Color','k')
            ylabel('SN','FontSize',12,'FontWeight','bold','Color','k')
            hold off
            hold on
            subplot(2,2,3)
            plot(rma(:,2),-rma(:,4),'k');
            xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
            ylabel('Depth UTM (WGS84)','FontSize',12,...
                'FontWeight','bold','Color','k')
            hold off
            
            %Function creating the parameters for the inversion matrix
            [lunpar, blocch, lunto, s, modv] =...
                segments_single(modv,origin(3),um,ma1);
            lunparz(1:length(lunpar),i)=lunpar;
            blocchi(1:length(blocch),i)=blocch;
            luntot(i)=lunto;
            sb(1:length(s),i)=s;
            
        end
    end
    
    pktime = tempi(i,PorS); %Picked time from SAC file
    
    if pa<3
        if isequal(f0,'[]')
            if evst==1
                ttheory(i,1)=D(i)/vth;
            elseif evst==2
                ttheory(i,1)=D(i)*degorutm/vth;
            end
            tempi(i,1)=pktime-ttheory(i,1);
        else
            tempi(i,1)=eval(f0);
        end   
    elseif pa==3
        ttheory(i,1)=sum(lunparz(:,i).*sb(:,i));
        tempi(i,1)=pktime-ttheory(i,1);
    end
       
    t00=tempis(1); %starting time of the waveform
    
    %starting-sample direct window
    cursor1 = floor((tempi(i,PorS)-t00)*srate);
    
    %Look for peak delay
    pdcursor=floor(cursor1+maxtpde*srate); %end sample for peak-delay pick
    tpdsisma = sisma(cursor1:pdcursor);% tapering direct wave
    fpdsisma = filter(Hd,tpdsisma);% filtering to find peak
    hpdsp = hilbert(fpdsisma); %hilbert
    pdsp = smooth(abs(hpdsp),nf/cf*srate);%ms of the filtered waveform
    [mspm,tspm]=max(pdsp);
    peakd(i,1)=tspm/srate;

    %Compute Qc from late lapse times (tCm)
    cursorc3 = floor((tempi(i,1)-t00+tCm-1)*srate);%coda start sample
    cursorc4 = floor(cursorc3 + tWm*srate-1);%coda end sample
    lsis=length(sisma);
    if cursorc4 > lsis
        tsismacm = sisma(cursorc3:lsis);%in case coda is too long
        indexlonger=indexlonger+1;
    else
        tsismacm = sisma(cursorc3:cursorc4);%tapering
    end
    L=length(tsismacm);
    tu=tukeywin(L,.05);
    tsismacm=tu.*tsismacm;    
    fsismacm = filter(Hd,tsismacm); %filter coda
    hspcm = hilbert(fsismacm); %hilbert
    spcm = smooth(abs(hspcm),nf/cf*srate);%ms of the filtered waveform
    lspm = length(spcm);
    tm = (tCm+1/srate:1/srate:tCm+lspm/srate)';
    
    %Only evaluate central time series
    edgeno=floor(0.05*length(tm));
    tm1=tm(edgeno:end-edgeno);
    spcm1=spcm(edgeno:end-edgeno);
    
    if nonlinear==0

        %THIS IS THE DATA VECTOR WITH LINEARISED THEORY - STANDARD
    
        EWz=spcm1.*tm1.^sped; 
        lspmz = log(EWz)/2/pi/cf; % source-station data
        Rz=corrcoef([tm1,lspmz]); %sets uncertainty
        polyz = polyfit(tm1,lspmz,1); %Qc^-1
    
        if polyz(1)<0
            Qm(i,1)=-polyz(1);
            RZZ(i,1)=abs(Rz(1,2));
        end    
    
    elseif nonlinear==1

    %THIS IS THE SYSTEM OF EQUATIONS FOR THE NON-LINEAR SOLUTION - USING
    %THE LINEAR INVERSE QC AS STARTING MODEL
        d_obs=zeros(ntW,1);
        ntm=0;
        for k=1:ntW
            lntm=length(ntm)+nW*srate;
%             if lntm>lspm
%                 break
%             end
            ntm = (k-1)*nW*srate + 1:k*nW*srate;
        
            nt= tm(floor(ntm));
            d_obs(k,1)=mean(spcm(floor(ntm)));    
        end
        d_obs1=d_obs(1:end-1)/d_obs(end);
    
        E=zeros(L1,1);    
        for n=1:L1
            d_pre=tlapse.^(-sped).*exp(-2*pi*cf.*tlapse*m1a(n));
            d_pre1=d_pre(1:end-1)/d_pre(end);
            E(n,1)=sum(abs(d_obs1-d_pre1));
        end
        [Emin, indexE] = min(E);
        Qm(i,1)=m1a(indexE);
        RZZ(i,1)=1/Emin;
    end
    
    if pa==3
        
        %Direct-wave intensity
        int = fin1*srate;% number of sample for the chosen window
        intn = fin1*srate;% number of sample for the noise
        cursor2 = floor(cursor1 + int-1); %end-sample of the direct window
        tsisma = sisma(cursor1:cursor2);% tapering direct wave
        fsisma = filter(Hd,tsisma);% filtering direct wave
        hsp = hilbert(fsisma); %hilbert
        sp = smooth(abs(hsp),nf/cf*srate);%ms of the filtered waveform
        spamp = trapz(sp); %direct energy
       
        %Coda measure for coda-normalization
        cursorc1 = floor((tempi(i,1)-t00+tCm-1)*srate);%coda start sample
        cursorc2 = floor(cursorc1 + tWm*srate-1);%coda end sample
        if cursorc2 > length(sisma)
            lsis=length(sisma);
            tsismac = sisma(cursorc1:lsis);%in case coda is too long
        else
            tsismac = sisma(cursorc1:cursorc2);%tapering
        end
        fsismac = filter(Hd,tsismac); %filter coda
        hspc = hilbert(fsismac); %hilbert
        spc = smooth(abs(hspc),nf/cf*srate);%ms of the filtered waveform
        spampc = trapz(spc); %coda energy
    
        %Noise
        if PorS == 2
            cursorn1 = floor((tempi(i,2)-t00-sn)*srate); %star sample noise
        elseif PorS == 3
            cursorn1 = floor((t00)*srate+1); %starting sample for noise
        end
        cursorn2 = floor(cursorn1 + intn-1); % end-sample noise-window
        tsisman = sisma(cursorn1:cursorn2);%tapering noise
        fsisman = filter(Hd,tsisman);%filtering noise
        hspn = hilbert(fsisman); %hilbert
        spn = smooth(abs(hspn),nf/cf*srate);%ms of the filtered waveform
        spampn = trapz(spn); %noise energy
    
        rmes = spamp/spampc;%spectral ratios signal-coda
        rmcn = spampc/spampn;%spectral ratios coda-noise
    
        %STORE
        signal(index,1) = spamp; %direct-wave energies
        coda(index,1) = spampc; %coda-wave energies
        rapsp(index,1) = rmes; %spectral ratios
        rapspcn(index,1)= rmcn; %coda-to-noise ratios
    
    end
end

noQm=sum(Qm==0)/length(Qm)*100;
if nonlinear==0
    noRZZ=sum(RZZ<=RZZ2)/length(RZZ)*100;
elseif nonlinear == 1
    noRZZ=sum(RZZ>=RZZ2)/length(RZZ)*100;
end
if noQm>20
    Qless=['Attention: ',num2str(noQm),...
        '% of your Q are <=0 and ',num2str(noRZZ),...
        '% of your uncertainty coefficients are lower than treshold'];
    
    warning(Qless)
    
else
    if noRZZ>20
        Rless=['Attention: ',num2str(noRZZ),...
            '% of your coerrelation coefficients are lower than treshold'];
        warning(Rless)
    end
end
    
if evst==2
    figure(rays)
    if degorutm==111 && hasMT==1
        load coastlines
        hold on
        geoshow(coastlat,coastlon);
        xlim([origin(1) origin(1)+nxc*stepg]);
        ylim([origin(2) origin(2)+nyc*stepg]);
    end
    if pa<3
        
        for nn=1:lls
            hold on
            plot([evestaz(nn,1) evestaz(nn,3)],...
                [evestaz(nn,2) evestaz(nn,4)],'k-')
        end
        hold on
        scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
        hold on
        scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
        hold off
        grid on
        ax = gca;
        ax.GridLineStyle = '-';
        ax.GridColor = 'k';
        ax.GridAlpha = 1;
        ax.LineWidth = 1;
        
        Murat.inversion.APeakDelay=Apd;
        Murat.inversion.AQCoda=Ac;
        Murat.inversion.epicentralDistance=D;
        Murat.geometry.even = even;
        Murat.geometry.station = staz;
        
    elseif pa==3
        
        Murat.geometry.even = even;
        Murat.geometry.station = staz;
        
        hold on
        
        subplot(2,2,1)        
        scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold on
        scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold off
        hold on
        subplot(2,2,2)        
        scatter(-even(:,3),even(:,2),sz,'c','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold on
        scatter(-staz(:,3),staz(:,2),sz,'^','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold off
        hold on
        subplot(2,2,3)        
        scatter(even(:,1),even(:,3),sz,'c','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold on
        scatter(staz(:,1),staz(:,3),sz,'^','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold off
    
    %INVERSION MATRIX for direct waves: A - in case of event-station in SAC
    % Only store the blocks crossed by at least 1 ray
        fover = modv(:,5)>=1; 
        inblocchi = modv(fover,:);
        inblocchi(:,5)=find(fover);

        A = zeros(length(Ac(:,1)),length(inblocchi(:,1)));
        for j = 1:length(sb(1,:))
            for i = 2:length(sb(:,1))
                bl = blocchi(i,j);
                slo = sb(i,j);
                l = lunparz(i,j);
                fbl = find(inblocchi(:,5)==bl);
                if isempty(fbl)==0
                    A(j,fbl)=-l*slo;% inversion matrix element in sec
                end
            end    
        end
    
        Murat.inversion.APeakDelay=Apd;
        Murat.inversion.AQCoda=Ac;
        Murat.inversion.epicentralDistance=D;
        Murat.inversion.AQ=A;
        
        Murat.inversion.partialLengths=lunparz;
        Murat.inversion.totalLengths=luntot;
        Murat.inversion.hitBlocks=inblocchi;
        Murat.inversion.slownessBlocks=sb;
        Murat.inversion.hitBlocks=inblocchi;
        
        if createrays == 1          
            Murat.results.rays=ma1;
        end
    
    end
    
    FName = 'Rays';
    saveas(rays,fullfile(FPath, FLabel, FName), fformat);
    
    if pa>1
        
        Qcsen=figure('Name','Qc sensitivity, first source-station pair',...
            'NumberTitle','off','visible',visib);
        Qcss=Ac(30,:);
        Qcs=zeros(nxc,nyc);

        index=0;
        for i=1:length(x)
            for j=1:length(y)
                index=index+1;        
                Qcs(i,j)=Qcss(index);
            end
        end
        
        [X,Y]=meshgrid(x,y);
        contourf(X,Y,Qcs')
        
        colorbar
        grid on
        ax = gca;
        ax.GridLineStyle = '-';
        ax.GridColor = 'k';
        ax.GridAlpha = 1;
        ax.LineWidth = 1;
        if degorutm==111 && hasMT==1
            hold on
            geoshow(coastlat,coastlon);
            xlim([origin(1) origin(1)+nxc*stepg]);
            ylim([origin(2) origin(2)+nyc*stepg]);
        end
        
        FName = 'Qc_sensitivity';
        saveas(Qcsen,fullfile(FPath, FLabel, FName), fformat);
    end
    
end  

% Setting up the inversion vector in case of 2- and 3-components data
icomp = 0;
ll=lls;

%Matrix A
if evst==2 && compon>1
    lA=length(A(:,1));
    Ac=Ac(1:compon:lA-compon,:);
    Apd=Apd(1:compon:lA-compon,:);
    if pa==3
        A=A(1:compon:lA-compon,:);
    end
    Murat.inversion.APeakDelay=Apd;
    Murat.inversion.AQCoda=Ac;
    Murat.inversion.AQ=A;
end

% 2 components (tipically the two horizontals)
if compon ==  2
    lsig = ll/2;
    signal1 = zeros(lsig,1); % The average direct wave energies
    coda1 = zeros(lsig,1); % The average coda wave energies
    rapsp1 = zeros(lsig,1); % The average spectral ratios
    rapspcn1 = zeros(lsig,1); % The average coda versus noise ratios
    Qm1 = zeros(lsig,1); % The average coda versus noise ratios
    RZZ1 = zeros(lsig,1); % The average coda versus noise ratios
    peakd1 = zeros(lsig,1); % The average coda versus noise ratios
    for i = 1:2:(ll-1)
        icomp = icomp+1;
        signal1(icomp,1) = (signal(i,1)+signal(i+1))/2;
        coda1(icomp,1) = (coda(i)+coda(i+1))/2;
        rapsp1(icomp,1) = (rapsp(i)+rapsp(i+1))/2;
        rapspcn1(icomp,1) = (rapspcn(i)+rapspcn(i+1))/2;
        Qm1=(Qm(i,1)+Qm(i+1))/2;
        RZZ1=(Rzz(i,1)+Rzz(i+1))/2;
        peakd1=(peakd(i,1)+lpeakd(i+1))/2;
    end    
    signal=signal1;
    coda=coda1;
    rapsp=rapsp1;
    rapspcn=rapspcn1;
    Qm=Qm1;
    RZZ=RZZ1;
    peakd=peakd1;
    
% 3 components (WE, SN, Z)
elseif compon == 3
    lsig = ll/3;
    signal1 = zeros(lsig,1); % The average direct wave energies
    coda1 = zeros(lsig,1); % The average coda wave energies
    rapsp1 = zeros(lsig,1); % The average spectral ratios
    rapspcn1 = zeros(lsig,1); % The average coda versus noise ratios
    Qm1 = zeros(lsig,1); % 
    RZZ1 = zeros(lsig,1); % 
    peakd1 = zeros(lsig,1); %
    
    for i = 1:3:(ll-2)
        icomp = icomp+1;
        signal1(icomp,1) = ((signal(i)+signal(i+1))/2 + signal(i+2))/2;
        coda1(icomp,1) = ((coda(i)+coda(i+1))/2 + coda(i+2))/2;
        rapsp1(icomp,1) = ((rapsp(i)+rapsp(i+1))/2 + rapsp(i+2))/2;
        rapspcn1(icomp,1) = ((rapspcn(i)+rapspcn(i+1))/2 + rapspcn(i+2))/2;
        Qm1(icomp,1) = ((Qm1(i)+Qm1(i+1))/2 + Qm1(i+2))/2;
        RZZ1(icomp,1) = ((RZZ1(i)+RZZ1(i+1))/2 + RZZ1(i+2))/2;
        peakd1(icomp,1) = ((peakd(i)+peakd(i+1))/2 + peakd(i+2))/2;
    end
    signal=signal1;
    coda=coda1;
    rapsp=rapsp1;
    rapspcn=rapspcn1;
    Qm=Qm1;
    RZZ=RZZ1;
    peakd=peakd1;
end

if pa<3
    lu=D;
    if evst==1
        time0=D/vth;
    elseif evst==2
        time0=D*degorutm/vth;
    end
    
elseif pa==3
    %Calculate average geometrical spreading factor and Q
    %sets the constant in the CN method
    lu=luntot;
    
    constCNQc=(tCm+tWm/2)^sped.*exp(Qm*2*pi*cf*(tCm+tWm/2));
    mcn = find(rapspcn<tresholdnoise);% set the weigth

    % weighting
    time0=tempi(1:compon:end-compon+1,PorS)-tempi(1:compon:end-compon+1,1);
    W1=rapspcn<tresholdnoise;
    W2=rapsp>1000;
    W3=rapsp<tresholdnoise;
    W4=W3+W2+W1;
    W=diag(1./(W4+1));% weights
    d0=log(rapsp./constCNQc)/2/pi/cf; %data of the inverse problem

    dW=W*d0;
    G1=-log(lu)/pi/cf; %matrix creation
    G1(:,2)=-time0;

    G=W*G1;
    %The 3 parameters are the constant, the geometrical spreading, and the
    %average Q and they are contained in constQmean.
    constQmean(1,1:2)=lsqlin(G,dW(:,1));% damped least square inversion
    cova = (G'*G)^(-1)*G'*cov(dW)*G*(G'*G)^(-1); %covariance matrix
    er = sqrt(diag(cova)); %error from the covariance matrix
    constQmean(2,:)=er;

end

%Peak-delay, for reference see e.g. Calvet et al. 2013, Tectonophysics
nop=find(peakd<1/cf);
peakd(nop)=[];
lu(nop)=[];
%Ac(nop,:)=[];

l10l=log10(lu);
l10pd=log10(peakd);
fitrobust = fit(l10l,l10pd,'poly1','Robust','on');
pab=[fitrobust.p1 fitrobust.p2];
l10=linspace(min(l10l),max(l10l));
lpd=pab(1)*l10+pab(2);
l10pdt=polyval(pab,l10l);
lpdelta=l10pd-l10pdt;

% To remove outliers
I = abs(lpdelta) > 2*std(lpdelta) | peakd >= maxtpde | peakd < mintpde;
outlierspd = excludedata(l10l,lpdelta,'indices',I);

%Same for Qc
%Remove anomalous Qm remembering they have a log-normal distribution
if nonlinear==0
    retainQm=find(Qm>0 & RZZ>RZZ2);
    mQm=mean(Qm(retainQm));
    outliersc = Qm>0 & Qm > mQm+2*std(Qm(retainQm));
    retainQm=find(Qm>0 & RZZ>RZZ2 & outliersc==0);
    discardQm=find(Qm<=0 | RZZ<RZZ2 | outliersc==1);
    discardQm2= find(Qm>0 & (RZZ<RZZ2 | outliersc==1));
    
elseif nonlinear==1
    retainQm=find(Qm>0 & RZZ<RZZ2);
    mQm=mean(Qm(retainQm));
    outliersc = Qm>0 & Qm > mQm+2*std(Qm(retainQm));
    retainQm=find(Qm>0 & RZZ<RZZ2 & outliersc==0);
    discardQm=find(Qm<=0 | RZZ>RZZ2 | outliersc==1);
    discardQm2= find(Qm>0 & (RZZ>RZZ2 | outliersc==1));
end

Qmdis=Qm;
Qm(discardQm)=mQm;

%save quantities for the data vector
Murat.data.measuredQc=Qm;
Murat.data.uncertainty=RZZ;
Murat.data.logPeakDelay=lpdelta;
Murat.data.outliersPeakDelay=outlierspd;
Murat.data.theoreticalTravelTime=time0;
Murat.data.averageQc=mQm;

if pa==3
    Murat.data.bodyEnergy=signal;
    Murat.data.codaEnergy=coda;
    Murat.data.ratioEnergy=rapsp;
    Murat.data.noiseCodaEnergy=rapspcn;
    Murat.data.averageQ=constQmean;
    Murat.data.dataQ=d0;
    Murat.data.weightQ=W;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS 

%plot to check that Qc is constant with travel time and peak delays
%increase with travel time

Qcpd=figure('Name','Qc and peak-delays',...
    'NumberTitle','off','visible',visib);
subplot(2,1,1)
plot(time0(retainQm),Qm(retainQm),'o',...
    'MarkerSize',6,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
hold on
plot(time0(discardQm2),Qmdis(discardQm2),'o',...
    'MarkerSize',6,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0])
hold off
title('Dependence of inverse Qc on travel time');
xlabel('Travel time (s)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Inverse Qc','FontSize',12,...
    'FontWeight','bold','Color','k')

subplot(2,1,2)
plot(fitrobust,'k--',l10l,l10pd,'ko',outlierspd,'r*')
title('Dependence of peak delays on travel times');
xlabel('Log. Travel time (s)','FontSize',12,'FontWeight','bold',...
    'Color','k')
ylabel('Log. Peak delay (s)','FontSize',12,...
    'FontWeight','bold','Color','k')
FName = 'Qc_Peak_Delay';
saveas(Qcpd, fullfile(FPath, FLabel, FName), fformat);

if pa==3
% The cyan dots represent the left-hand side of the
% CN method linear equation. The red line is the fit of the three
% parameters. The black lines are the uncertainties.
% In the lower panel the coda-to noise ratios, to see the effe ct of noise
% on the data.

    % Residual geometrical spreading in the coda
    gc=polyfit(time0,log(rapspcn),1);
    totgs=-(constQmean(1,1));
    totgs1=-(constQmean(1,1)-er(1));
    totgs2=-(constQmean(1,1)+er(1));
    % Right hand side of the CN linear equation - corrected by residual
    % geometrical spreading
    x = totgs*log(luntot)/pi/cf...
        -time0*constQmean(1,2);
    xer1 = totgs1*log(luntot)/pi/cf...
        -time0*(constQmean(1,2)-er(2));
    xer2 = totgs2*log(luntot)/pi/cf...
        -time0*(constQmean(1,2)+er(2));

%Plot of the left and right hand sides of the CN equation. The average
%inverse Q and geometrical spreading coefficient are schown above each
%panel
    CN=figure('Name','Coda-normalized intensities vs travel times',...
        'NumberTitle','off','visible',visib);
    subplot(2,1,1)
    plot(time0,d0, 'o',...
        'MarkerSize',6,'MarkerEdgeColor',[0 0.6 1],...
        'MarkerFaceColor',[0 0.6 1])
    hold on
    plot(time0,x,'r.')
    plot(time0,xer1,'k.')
    plot(time0,xer2,'k.')
    hold off
    xlabel('Travel time (s)','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('Logarithm of the direct-to-coda energy ratio','FontSize',12,...
        'FontWeight','bold','Color','k')

    title(['Average inverse quality factor = ',...
        num2str(constQmean(1,2)),' +/- ', num2str(er(2)),...
        ' and geometrical spreading = ',...
        num2str(constQmean(1,1)),' +/- ', num2str(constQmean(2,1))],...
        'FontSize',12,'FontWeight','bold','Color','k')
    h_legend=legend('Direct-to-coda energy ratios','Inverse average Q',...
        'Inverse Q spreading-related uncertainties');
    set(gca,'XTick',0:round(max(time0))) ;
    set(gca,'YTick',min(log(rapsp/constCNQc)+...
        constQmean(1,1)*log(luntot))/pi/cf:...
        (max(log(rapsp/constCNQc)+...
        constQmean(1,1)*log(luntot))/pi/cf-...
        min(log(rapsp/constCNQc)+constQmean(1,1)*log(luntot))/pi/cf)/5:...
        max(log(rapsp/constCNQc)+constQmean(1,1)*log(luntot))/pi/cf) ;
    set(gca,'FontSize',10,'FontWeight','bold');
    set(h_legend,'FontSize',10,'FontWeight','bold');

% %Plot of coda-to-noise energy ratios. The residual geometrical
% spreading should be around zero, or coherent phases could be
% included in the coda-energy; it also provides the average coda-to-noise
% ratio
    subplot(2,1,2)
    plot(time0,(log(rapspcn)),'o','MarkerSize',6,'MarkerEdgeColor',...
        [.2 .8 0],'MarkerFaceColor',[.2 .8 0])
    mc = mean(rapspcn);
    mno=length(mcn);
    perno= mno/length(rapspcn)*100;
    hold on
    plot(time0(mcn),log(rapspcn(mcn)), 'o','MarkerSize',6,...
    'MarkerEdgeColor',[1 .5 0],'MarkerFaceColor',[1 .5 0])
    hold off
    xlabel('Travel time (s)','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('Logarithm of the coda-to-noise energy ratio','FontSize',12,...
        'FontWeight','bold','Color','k')
    title(['Residual coda geometrical spreading = ',num2str(gc(1)), ...
        ' while the percentage of measures < ',num2str(tresholdnoise),...
        ' is ',num2str(perno), '%'],'FontSize',12,...
        'FontWeight','bold','Color','k');
    h_legend1=legend('Coda-to-noise energy ratios',...
        'Ratios below treshold');
    set(gca,'XTick',0:round(max(time0)));
    set(gca,'YTick',min(log(rapspcn)):...
        (max(log(rapspcn))-min(log(rapspcn)))/5:...
        max(log(rapspcn)));
    set(gca,'FontSize',10,'FontWeight','bold');
    set(h_legend1,'FontSize',10,'FontWeight','bold');

    text(max(time0), min(log(rapspcn)),...
    [' Lapse time = ',num2str(tCm),' s, coda window = ',...
    num2str(tWm), ' s'],'VerticalAlignment','bottom',...
    'HorizontalAlignment','right','FontSize',12)
    FName = 'CN_plot';
    saveas(CN,fullfile(FPath, FLabel, FName), fformat);
end
% END PLOTS

clearvars -except Murat
save('Murat.mat','Murat');

% rtoc=toc;
%%  2D peak-delay and Qc, and 3D CN TOMOGRAPHIC INVERSIONS
clear
load Murat.mat
pa=Murat.analysis;

%PATHS and FIGURES
FPath=Murat.paths.workingdir;
FLabel=Murat.paths.label;
fformat=Murat.figures.format;

%DATA
lls=length(Murat.paths.listasac);
cf=Murat.data.centralFrequency;
sped=Murat.data.spectralDecay;
Qm=Murat.data.measuredQc;
RZZ=Murat.data.uncertainty;
lpdelta=Murat.data.logPeakDelay;
outlierspd=Murat.data.outliersPeakDelay;
time0=Murat.data.theoreticalTravelTime;
mQm=Murat.data.averageQc;

%GEOMETRY
evst=Murat.geometry.import;
XY=Murat.geometry.map;
nxc=Murat.geometry.gridX;
nyc=Murat.geometry.gridY;


% INVERSION
sizea=Murat.inversion.sizeCheck;
latt=Murat.inversion.lowCheck;
hatt=Murat.inversion.highCheck;
Apd=Murat.inversion.APeakDelay;
Ac=Murat.inversion.AQCoda;
D=Murat.inversion.epicentralDistance;
        
if pa==3
    modv=Murat.geometry.modv;
    resol2=Murat.geometry.resolutionMin;
    dtreshold=Murat.geometry.depthTreshold;
    signal=Murat.data.bodyEnergy;
    coda=Murat.data.codaEnergy;
    rapsp=Murat.data.ratioEnergy;
    rapspcn=Murat.data.noiseCodaEnergy;
    constQmean=Murat.data.averageQ;
    d0=Murat.data.dataQ;
    W=Murat.data.weightQ;
    if evst==1
        A=Murat.inversion.AQ;
        lunparz=Murat.inversion.partialLengths;
        luntot=Murat.inversion.totalLengths;
        inblocchi=Murat.inversion.hitBlocks;
        sb=Murat.inversion.slownessBlocks;
    end
end

pd=XY(:,1);
pd(:,2)=XY(:,2);
pd(:,3)=-1000;
Qc=pd;

%Peak delay mapping
lpdelta_o=lpdelta(outlierspd==0);
Apd=Apd(outlierspd==0,:);
Apd(Apd~=0)=1;
lApd=size(Apd);
Apd1=Apd;
pdelay=zeros(lApd(2),1);
for i=1:lApd(1)
    Apd1(i,1:end)=Apd1(i,1:end)*lpdelta_o(i);
end
for j=1:lApd(2)
    a1=Apd1(:,j);
    sipd=find(a1);
    pdelay(j,1)=mean(a1(sipd));
end
pdelay(isnan(pdelay))=0;
pd(:,4)=pdelay;

%Qc mapping
%Remove anomalous Qm
retainQm=find(Qm~=mQm);
stat=Qm(retainQm);

if pa==1

    p = [0.1 0.9];
    y = quantile(stat,p);
    zk = y;
    mm = [mQm median(stat)];
    sk = [skewness(stat) kurtosis(stat)];

    lAc=size(Ac);
    Ac1=Ac;

    Qcf=zeros(lAc(2),1);
    for i=1:lAc(1)
        Ac1(i,1:end)=Ac1(i,1:end)*Qm(i);
    end
    for j=1:lAc(2)
        a1=Ac1(:,j);
        sic=find(a1);
        Qcf(j,1)=mean(a1(sic));
    end
    Qcf(isnan(Qcf))=0;
    Qcf(Qcf==0)=mQm;
    Qc(:,4)=Qcf;

elseif pa>1
    Ac1=Ac(retainQm,:);
    RZZ1=RZZ(retainQm,1);
    W1=RZZ1<0.3;
    W2=RZZ1<0.5;
    W3=RZZ1<0.7;
    W4=W3+W2+W1;
    Wc=diag(1./(W4+1));% weights
    
    dcW=Wc*stat;
    Gc=Wc*Ac1;
    [lA,llA]=size(Gc);
    [Uc,Sc,Vc]=svd(Gc);

    index1=0;
    
    % smooting parameter is user defined
    LcQc=figure('Name','L-curve Qc','NumberTitle','off');
    l_curve(Uc,diag(Sc),stat,'Tikh')
    tik0_regC=input('Your personal smoothing parameter for coda ');
    FName = 'Lc_Qc';
    saveas(LcQc,fullfile(FPath, FLabel, FName), fformat);
    close(LcQc)
    
    % picard plot
    PpQc=figure('Name','Picard-plot Qc',...
        'NumberTitle','off','visible','off');
    picard(Uc,diag(Sc),stat);
    FName = 'Picard_Qc';
    saveas(PpQc,fullfile(FPath, FLabel, FName), fformat);
    
    mtik0C=tikhonov(Uc,diag(Sc),Vc,stat,tik0_regC);
    Qc(:,4)=mtik0C;

    %Testing - Creating 2D checkerboard matrix
    nxc1=nxc/sizea;
    nyc1=nyc/sizea;
    I = imresize(xor(mod(1:nyc1, 2).', mod(1:nxc1, 2)),...
        [sizea*nyc1 sizea*nxc1], 'nearest');
    Qc(:,5)=I(1:end);
    Qc(Qc(:,5)==1,5)=latt;
    Qc(Qc(:,5)==0,5)=hatt;
    
    Qc5= Qc(:,5);
    re = Gc*Qc5;
    mcheckc=tikhonov(Uc,diag(Sc),Vc,re,tik0_regC);
    Qc(:,6)=mcheckc;
    Qc(:,7)=Gc(1,:);
    Qc(:,8)=Gc(2,:);

end

if pa==3
% Direct wave attenuation

    %Data creation, removing the pre-calculated parameters
    d1 = d0  + constQmean(1,1)*log(luntot)/pi/cf...
    + time0*constQmean(1,2);


    % if the ray crosses a block not solved by the inversion, the
    % block is characterized by the average quality factor, and the data
    % vector is updated.
    
    A1=A;
    A=W*A1;

    % NEW DATA VECTOR

    dW1=W*d1;

    % tikhonov inversion - by using the programs in HANSEN et al. 1994
    [U,S,V]=svd(A);

    %sets the smoothing parameter - always user defined
    LcCN=figure('Name','L-curve coda-normalization','NumberTitle','off');
    l_curve(U,diag(S),dW1,'Tikh')
    tik0_reg=input('Your personal smoothing parameter ');
    FName = 'Lc_CN';
    saveas(LcCN,fullfile(FPath, FLabel, FName), fformat);
    close(LcCN)
    
    % picard plot
    PpCN=figure('Name','Picard coda-norm.','NumberTitle','off',...
        'visible','off');
    picard(U,diag(S),dW1);
    FName = 'Picard_CN';
    saveas(PpCN,fullfile(FPath, FLabel, FName), fformat);
    
    %results
    mtik0=tikhonov(U,diag(S),V,dW1,tik0_reg);
    mQ1=[inblocchi(:,1:3) mtik0];

    % Simple resolution tests: produces columns in output.
    %FIRST TEST: CHECKERBOARD
    %INPUT: Checkerboard structure - node spacing of the anomalies
    % doubled with respect to the grid.

    %Depth
    passox = find(modv(:,1)~=modv(1,1),1,'first')-1;
    passoy = find(modv(:,2)~=modv(1,2),1,'first')-1;

    sizea2=2*sizea;
    sizeap=sizea*passoy;
    sizeap1=(sizea+1)*passoy;

    for k=1:sizea
        modv(k:sizea2:passoy-sizea+k,6)=hatt;
        modv(sizea+k:sizea2:passoy-sizea+k,6)=latt;
    end
    for k=1:sizea-1
        modv(k*passoy+1:(k+1)*passoy,6)=modv(1:passoy,6);
    end
    for k=1:sizea
        modv(sizeap+k:sizea2:sizeap1-sizea+k,6)=latt;
        modv(sizeap+sizea+k:sizea2:sizeap1-sizea+k,6)=hatt;
    end

    py4  = 2*sizeap;
    for k=1:sizea-1
        modv((sizea+k)*passoy+1:(sizea+k+1)*passoy,6)=...
            modv(sizeap+1:sizeap1,6);
    end
    z = (passox-mod(passox,py4))/py4;
    for i = 1:(z-1)
        modv(i*py4+1:(i+1)*py4,6)=modv(1:py4,6);
    end
    if ~isequal(mod(passox,py4),0)
        modv(z*py4+1:z*py4+mod(passox,py4),6)= modv(1:mod(passox,py4),6);
    end

    %Along y
    sizeapx=sizea*passox;

    for k=1:sizea-1
        modv(k*passox+1:(k+1)*passox,6)=modv(1:passox,6);
    end

    for k = 1:sizeapx
        if modv(k,6)==hatt
            modv(sizeapx+k,6)=latt;
        elseif modv(k,6)==latt
            modv(sizeapx+k,6)=hatt;
        end
    end

    %Along x
    px4  = 2*sizea*passox;

    z2= (length(modv(:,1))-mod(length(modv(:,1)),px4))/px4;
    for i = 1:(z2-1)
        modv(i*px4+1:(i+1)*px4,6)=modv(1:px4,6);
    end
    if ~isequal(mod(passox,py4),0)
        modv(z*px4+1:z*px4+mod(length(modv(:,1)),px4),6)=...
            modv(1:mod(length(modv(:,1)),px4),6);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RESULTS OF THE FIRST TOMOGRAPHIC INVERSION + tests
% In the original UTM reference system we put the quality factors in
% the center of each block. The average quality factor characterizes the
% blocks not solved in the inversion.
    Q3D=modv(:,1:3);
    Q3D(:,4)=constQmean(1,2);
    % create a file of input
    in1=zeros(length(mtik0),1);
    input01 = zeros(length(mtik0),1);

    for i = 1:length(in1)
        in1(i)=find(modv(:,1)==inblocchi(i,1)...
            & modv(:,2)==inblocchi(i,2)...
            & modv(:,3)==inblocchi(i,3));
        Q3D(in1(i),4)=mtik0(i)+constQmean(1,2);
        input01(i)=modv(in1(i),6);
    end

    % result in the original reference system
    Q3D(:,1)=modv(:,1)+resol2;
    Q3D(:,2)=modv(:,2)+resol2;
    Q3D(:,3)=modv(:,3)+resol2;

    % input save
    Qin=[Q3D(:,1:3) modv(:,6)];

    A1=W*A;
    % output file and save
    re = A1*input01;
    Qou(:,1:3)=Q3D(:,1:3);
    Qou(:,4)=0;
    mcheck=tikhonov(U,diag(S),V,W*(re-mean(input01)),tik0_reg);
    Qou(in1,4)=mcheck;

    % WE ALSO OUTPUT THE DIAGONAL OF THE RESOLUTION MATRIX
    % using the filter functions
    sS=size(S);
    fil_reg = fil_fac(diag(S),tik0_reg);
    if sS(2)>sS(1)
        fil_reg(end+1:sS(2),1)=0;
    end
    R=V*diag(fil_reg)*V';
    dR =diag(R);
    Qou(in1,5)=dR;

    % AS THIRD AND FINAL TEST, THE USER CAN SET SYNTHETIC ANOMALIES.
    %First is with the entire set of actual results - this is improper but
    %used in some studies
    lmQ1=length(mQ1(:,4));
    syQ0=constQmean(1,1)*log(luntot)/pi/cf+time0.*constQmean(1,2)+...
        A*mQ1(:,4);

    G1=-log(luntot)/pi/cf; %matrix creation
    G1(:,2)=-time0;

    dsW=W*syQ0;
    G=W*G1;

    %The 2 parameters are the geometrical spreading and the
    %average Q and they are contained in synthQmean.
    synthQmean(1,1:2)=lsqlin(G,dsW(:,1));% damped least square inversion
    cova = (G'*G)^(-1)*G'*cov(dsW)*G*(G'*G)^(-1); %covariance matrix
    er = sqrt(diag(cova)); %error from the covariance matrix
    synthQmean(2,:)=er;

    %data creation
    syQ1 = W*(syQ0 + synthQmean(1,1)*log(luntot)/pi/cf...
        + time0*synthQmean(1,2));

    [syU,syS,syV]=svd(A);
    
    %results
    sytik0=tikhonov(syU,diag(syS),syV,syQ1,tik0_reg);
    Qou(:,6)=synthQmean(1,2);
    Qou(in1,6)=synthQmean(1,2)+sytik0;
    Q3D(:,5)=Qin(:,4);
    Q3D(:,6:8)=Qou(:,4:6);

    %Second is user-defined - better
    %As example, create a 2-layer medium with average Q latt and hatt
    Qplus=find(inblocchi(:,3)>dtreshold);
    Qminus=find(inblocchi(:,3)<dtreshold);
    A2=A;
    syQ2=mQ1;
    syQ2(Qplus,4)=latt;
    syQ2(Qminus,4)=hatt;
    A3=zeros(lls,1);
    for k=1:lls
        A3(k,1)=A2(k,:)*syQ2(:,4);
    end
        

    tsp=sped; %theoretical geom. spr.
    d0=tsp*log(luntot)/pi/cf-A3;
    dW=W*d0;
    G1=-log(luntot)/pi/cf;
    G1(:,2)=-time0;

    G=W*G1;
    syconstQmean(1,1:2)=lsqlin(G,dW(:,1));% damped least square inversion
    cova = (G'*G)^(-1)*G'*cov(dW)*G*(G'*G)^(-1); %covariance matrix
    er = sqrt(diag(cova)); %error from the covariance matrix
    syconstQmean(2,:)=er;

    syQ3=tsp*log(luntot)/pi/cf-A3+time0*constQmean(1,2);

    dsW2=W*syQ3;
    
    sytik02=tikhonov(syU,diag(syS),syV,dsW2,tik0_reg);

    inup=Qou(:,3)>dtreshold;
    indown=Qou(:,3)<dtreshold;
    Qou(inup,7)=hatt;
    Qou(indown,8)=latt;
    Qou(in1,7)=syQ2(:,4);
    Qou(in1,8)=syconstQmean(1,2)+sytik02;
    Q3D(:,9:10)=Qou(:,7:8);

    % save Q
    FName = 'Q3D.txt';
    save(fullfile(FPath, FLabel, FName), 'Q3D','-ascii');
end
% save peak-delay
FName = 'peakdelay.txt';
save(fullfile(FPath, FLabel, FName), 'pd','-ascii');
% save Qc
FName = 'Qc.txt';
save(fullfile(FPath, FLabel, FName), 'Qc','-ascii');

%% Creating maps

clear

load Murat.mat
hasMT=Murat.figures.hasMT;

if hasMT
    load coastlines
end

pa=Murat.analysis;

%PATHS and FIGURES
FPath=Murat.paths.workingdir;
FLabel=Murat.paths.label;
visib=Murat.figures.visibility;
sz=Murat.figures.sizeMarker;
fformat=Murat.figures.format;

%GEOMETRY
origin=Murat.geometry.origin;
stepg=Murat.geometry.gridStep;
nxc=Murat.geometry.gridX;
nyc=Murat.geometry.gridY;
x=Murat.geometry.x;
y=Murat.geometry.y;
even=Murat.geometry.even;
staz=Murat.geometry.station;
degorutm=Murat.geometry.degreesorutm;

if pa==3
    modv=Murat.geometry.modv;
    um=Murat.geometry.unity;
    sections=Murat.figures.sections;
    % load Q3D
    FName_Q3D = 'Q3D.txt';
    Q3D=load(fullfile(FPath, FLabel, FName_Q3D));
end

% load peak-delay
FName_pd = 'peakdelay.txt';
pd=load(fullfile(FPath, FLabel, FName_pd));

% load Qc
FName_Qc = 'Qc.txt';
Qc=load(fullfile(FPath, FLabel, FName_Qc));

pdel=zeros(length(x),length(y));
QQc=zeros(length(x),length(y));
QQchi=zeros(length(x),length(y));
QQcho=zeros(length(x),length(y));
[X,Y]=meshgrid(x,y);
index=0;
for i=1:length(x)
    for j=1:length(y)
        index=index+1;        
        pdel(i,j)=pd(index,4);
        QQc(i,j)=Qc(index,4);
        if pa>1
            QQchi(i,j)=Qc(index,5);
            QQcho(i,j)=Qc(index,6);
        end
    end
end

pdmap=figure('Name','Peak-delay map','NumberTitle','off','visible',visib);
contourf(X,Y,pdel');
axis equal
view(2)
colormap(autumn);
hcb=colorbar;
title(hcb,'Log. Peak Delay','FontSize',14,'FontWeight','bold','Color','k');
    
xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
ylabel('SN','FontSize',12,'FontWeight','bold','Color','k')
title('Peak-delay variations',...
    'FontSize',12,'FontWeight','bold','Color','k');
hold on
scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
hold on
scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
if degorutm==111 && hasMT==1
    hold on
    geoshow(coastlat,coastlon);
    xlim([origin(1) origin(1)+nxc*stepg]);
    ylim([origin(2) origin(2)+nyc*stepg]);
end
hold off
FName = 'Peak_delay_map';
saveas(pdmap,fullfile(FPath, FLabel, FName), fformat);

Qcmap=figure('Name','Qc map','NumberTitle','off','visible',visib);
contourf(X,Y,QQc');
axis equal
view(2)
colormap(copper)

hcb=colorbar;
title(hcb,'Inverse Qc','FontSize',14,'FontWeight','bold','Color','k');
    
xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
ylabel('SN','FontSize',12,'FontWeight','bold','Color','k')
title('Coda attenuation variations','FontSize',12,'FontWeight','bold',...
    'Color','k');
hold on
scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
hold on
scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
if degorutm==111 && hasMT==1
    hold on
    geoshow(coastlat,coastlon);
    xlim([origin(1) origin(1)+nxc*stepg]);
    ylim([origin(2) origin(2)+nyc*stepg]);
end
hold off
FName = 'Qc_map';
saveas(Qcmap,fullfile(FPath, FLabel, FName), fformat);

%Parameter analysis
if pa ==1
    pdd = pd(:,4)~=0 & Qc(:,4)~=mQm;

    pdef=pd(pdd,:);
    Qcef=Qc(pdd,:);
    pd(:,5)=0;

elseif pa>=2
    
    Qcchecki=figure('Name','Qc checkerboard test input','NumberTitle',...
        'off','visible',visib);
    contourf(X,Y,QQchi');
    axis equal
    view(2)
    colormap(gray)
    
    hcb=colorbar;
    title(hcb,'Inverse Qc','FontSize',14,'FontWeight','bold','Color','k');
    
    xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('SN','FontSize',12,'FontWeight','bold','Color','k')
    title('Qc checherboard input',...
        'FontSize',12,'FontWeight','bold','Color','k');
    
    hold on
    scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    hold on
    scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    if degorutm==111 && hasMT==1
        hold on
        geoshow(coastlat,coastlon);
        xlim([origin(1) origin(1)+nxc*stepg]);
        ylim([origin(2) origin(2)+nyc*stepg]);
    end
    hold off
    FName = 'Qc_checkerboard_input';
    saveas(Qcchecki,fullfile(FPath, FLabel, FName), fformat);
    
    Qcchecko=figure('Name','Qc checkerboard test output','NumberTitle',...
        'off','visible',visib);
    contourf(X,Y,QQcho');
    axis equal
    view(2)
    colormap(gray)
    
    hcb=colorbar;
    title(hcb,'Inverse Qc','FontSize',14,'FontWeight','bold','Color','k');
    
    xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('SN','FontSize',12,'FontWeight','bold','Color','k')
    title('Qc checherboard output',...
        'FontSize',12,'FontWeight','bold','Color','k');
    
    hold on
    scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1.5)
    hold on
    scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1.5)
    if degorutm==111 && hasMT==1
        hold on
        geoshow(coastlat,coastlon);
        xlim([origin(1) origin(1)+nxc*stepg]);
        ylim([origin(2) origin(2)+nyc*stepg]);
    end
    hold off
    
    FName = 'Qc_checkerboard_output';
    saveas(Qcchecko,fullfile(FPath, FLabel, FName), fformat);
    
    %Parameter analysis
    pdd =  pd(:,4)~=0;
    pdef=pd(pdd,:);
    Qcef=Qc(pdd,:);
    
end

pdef(:,4)=pdef(:,4)-mean(pdef(:,4));
Qcef(:,4)=Qcef(:,4)-mean(Qcef(:,4));
mipdm=min(pdef(:,4));
mapdm=max(pdef(:,4));
miQcm=min(Qcef(:,4));
maQcm=max(Qcef(:,4));
trepd=0.05*(mapdm-mipdm)/2;
treQc=0.05*(maQcm-miQcm)/2;

Qps=Qcef(:,4);
pdps=pdef(:,4);

param_plot=figure('Name','Parameter space separation',...
    'NumberTitle','off','visible',visib);
hax=axes;

par=pdef(:,1:2);
par(:,3)=pdef(:,3);

c=Qps<-treQc & pdps<-trepd;
par(c,4)=1;
scatter(Qps(c),pdps(c),65,'filled','MarkerFaceColor',[0 0.8 0])
hold on
line([0 0],[mipdm-trepd mapdm+trepd],'Color',[0 0 0],...
    'LineWidth',3)
hold on
line([miQcm-treQc maQcm+treQc],[0 0],'Color',[0 0 0],...
    'LineWidth',3)
hold on
c=Qps<-treQc & pdps>trepd;
par(c,4)=2;
scatter(Qps(c),pdps(c),65,'filled','MarkerFaceColor',[0 0.6 1])
hold on
c=Qps>treQc & pdps<-trepd;
par(c,4)=3;
scatter(Qps(c),pdps(c),65,'filled','MarkerFaceColor',[1 0.6 0])
hold on
c=Qps>treQc & pdps>trepd;
par(c,4)=4;
scatter(Qps(c),pdps(c),65,'filled','MarkerFaceColor',[1 0 0])
hold on
c=(Qps>-treQc & Qps<treQc) | (pdps>-trepd & pdps<trepd);
par(c,4)=0;
scatter(Qps(c),pdps(c),85,'filled','MarkerFaceColor',[0.7 0.7 0.7],...
    'MarkerEdgeColor',[1 1 1],'LineWidth',2)
hold off 
xlim([miQcm-treQc maQcm+treQc])
ylim([mipdm-trepd mapdm+trepd])
xlabel('Qc','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Log. peak delay','FontSize',12,'FontWeight','bold','Color','k')
title('Parameter space plot',...
    'FontSize',12,'FontWeight','bold','Color','k');
FName = 'Parameter_space_variations';
saveas(param_plot,fullfile(FPath, FLabel, FName), fformat);

para=pd(:,1:3);
for k=1:length(par(:,1))
    px=par(k,1);
    py=par(k,2);
    pp=par(k,4);
    pf= pd(:,1)==px & pd(:,2)==py;
    para(pf,4)=pp;
end

index=0;
param=zeros(size(QQc));
for i=1:length(x)
    for j=1:length(y)
        index=index+1;        
        param(i,j)=para(index,4)-5;
    end
end

mparam=figure('Name','Parameter separation map',...
    'NumberTitle','off','visible',visib);
surf(X,Y,param');
view(2)
un_X = unique(param);

fu=find(un_X==-1);
if isempty(fu)
    cmap = [0.7 0.7 0.7;  0 0.8 0; 0 0.6 1; 1 0.6 0];
    HTick={'Average','Ls La','Hs La','Ls Ha'};
else
    cmap = [0.7 0.7 0.7;  0 0.8 0; 0 0.6 1; 1 0.6 0; 1 0 0];
    HTick={'Average','Ls La','Hs La','Ls Ha','Hs Ha'};
end
colormap(cmap)

colorbar('Ticks',un_X,'TickLabels',HTick);
    
hold on
scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1.5)
hold on
scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1.5)
if degorutm==111 && hasMT==1
    hold on
    geoshow(coastlat,coastlon);
    xlim([origin(1) origin(1)+nxc*stepg]);
    ylim([origin(2) origin(2)+nyc*stepg]);
end
hold off
xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
ylabel('SN','FontSize',12,'FontWeight','bold','Color','k')
axis square
title('Parameter separation','FontSize',12,'FontWeight','bold',...
    'Color','k');
FName = 'Parameter_map';
saveas(mparam,fullfile(FPath, FLabel, FName), fformat);

%Image velocity model
if pa==3
    xD=unique(modv(:,1));
    yD=unique(modv(:,2));
    zD=sort(unique(modv(:,3)),'descend');
    
    ix=length(xD);%numer of x layers,given the step of the grid
    iy=length(yD);%numer of y layers,given the step of the grid
    iz=length(zD);%numer of depths layers,given the step of the grid
    
    V=zeros(iy,ix,iz);
    Q=zeros(iy,ix,iz);
    Qi=zeros(iy,ix,iz);
    Qo=zeros(iy,ix,iz);
    index=0;
    for i=1:ix
        for j=1:iy
            for k=1:iz
                index=index+1;
                V(j,i,k)=modv(index,4);
                Q(j,i,k)=Q3D(index,4);
                Qi(j,i,k)=Q3D(index,5);
                Qo(j,i,k)=Q3D(index,6);
            end
        end
    end
    
    [X,Y,Z]=meshgrid(xD,yD,zD);
    
    V_map=figure('Name','Velocity Model','NumberTitle','off',...
        'visible',visib);
    WEi=sections(1,:);
    SNi=sections(2,:);
    zi=sections(3,:);
    %Velocity model
    slice(X,Y,Z,V,WEi,SNi,zi,'spline');
    colormap(cool)
    hcb=colorbar;
    title(hcb,'V','FontSize',14,'FontWeight','bold','Color','k');
    axis equal
    view(2)
    hold on
    scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    hold on
    scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    hold off
    xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('SN','FontSize',12,'FontWeight','bold','Color','k')
    zlabel('Depth (m)','FontSize',12,'FontWeight','bold','Color','k')
    title('Velocity model',...
        'FontSize',12,'FontWeight','bold','Color','k');
    
    FName = 'V_model';
    savefig(V_map,fullfile(FPath, FLabel, FName));
    
    Q_map=figure('Name','Q Model','NumberTitle','off',...
        'visible',visib);
    
    %Attenuation model
    slice(X,Y,Z,Q,WEi,SNi,zi,'spline');
    colormap(autumn)
    hcb=colorbar;
    title(hcb,'Inverse Q','FontSize',14,'FontWeight','bold','Color','k');
    axis equal
    view(2)
    hold on
    scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    hold on
    scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    hold off
    xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('SN','FontSize',12,'FontWeight','bold','Color','k')
    zlabel('Depth (m)','FontSize',12,'FontWeight','bold','Color','k')
    title('Attenuation model',...
        'FontSize',12,'FontWeight','bold','Color','k');
    
    FName = 'Q_model';
    savefig(Q_map,fullfile(FPath, FLabel, FName));
    
    Qchecki =...
        figure('Name','Checkerboard 3D input','NumberTitle','off',...
        'visible',visib);
    
    %3D checkerboard input
    slice(X,Y,Z,Qi,WEi,SNi,zi,'spline');
    colormap(gray)
    hcb=colorbar;
    title(hcb,'Inverse Q','FontSize',14,'FontWeight','bold','Color','k');
    axis equal
    view(2)
    hold on
    scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    hold on
    scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    hold off
    xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('SN','FontSize',12,'FontWeight','bold','Color','k')
    zlabel('Depth (m)','FontSize',12,'FontWeight','bold','Color','k')
    title('3D checkerboard input',...
        'FontSize',12,'FontWeight','bold','Color','k');
    
    FName = '3DQ_check_input';
    savefig(Qchecki,fullfile(FPath, FLabel, FName));
    
    Qchecko =...
        figure('Name','Checkerboard 3D output','NumberTitle','off',...
        'visible',visib);
    
    %3D checkerboard output
    slice(X,Y,Z,Qo,WEi,SNi,zi,'spline');
    colormap(gray)
    hcb=colorbar;
    title(hcb,'Inverse Q','FontSize',14,'FontWeight','bold','Color','k');
    
    axis equal
    view(2)
    hold on
    scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    hold on
    scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    hold off
    xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('SN','FontSize',12,'FontWeight','bold','Color','k')
    zlabel('Depth (m)','FontSize',12,'FontWeight','bold','Color','k')
    title('3D checkerboard output',...
        'FontSize',12,'FontWeight','bold','Color','k');
    
    FName = '3DQ_check_output';
    savefig(Qchecko,fullfile(FPath, FLabel, FName));
end

