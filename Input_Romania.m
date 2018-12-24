%% INPUTS
addpath('/Users/lucadesiena/Documents/Utilities_Matlab')
clear
close all
clc

disp('Input Section')

% The following prompts are included TO BE  EDITED before running!

% Which analysis do you want to perform?
% Pick delay and Qc without kernels: pa=1
% Pick delay and Qc with kernels: pa=2
% Pick delay, kernel-Qc and P/S wave attenuation with the CN method: pa=3
pa=2;

% Folder containing data
DFolder = './sac_Romania/*.sac';

% Creates folders and paths to store results and figures
FLabel = 'Romania';
FPath = './';

if exist(FLabel,'dir')~=7
    mkdir(FLabel)
end

% Storing inversion parameters
inputinv = 'inversion.mat';
inputdata = 'data.mat';

% Output figure format 
fformat='jpeg';

% Figure visibility - pick 'on' or 'off'
visib='on';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs necessary for all the analyses
% Choose between P- (2) and S- (3) waves for peak delay
% RULE: First column in time files always contains origin of events
% RULE: Second column in time files always contains P-wave phase
% RULE: Third column in time files always contains S-wave phase
PorS = 2;

% Work with 1 vertical (1) or 2 horizontal (2) recordings, or with the
% three components(3)
compon = 1;

% Length of the window used to measure noise and direct wave intensity
fin1 = 1;%edit

% To see spectrogram (set=1) of the seismogram seisi: set seisi='number of
% seismogram to see spectrogram'
seesp=0;
seisi=100;

% Central frequency - set it according to your spectrograms
cf = 1;

% Maximum window to pick pick-delays in seconds
maxtpde = 25;

% Minimum peak delay considering scattering in the area and frequency
mintpde = 0.5;%1/cf;

% Lapse time for the start of the window used to measure and calculate the
% normalization energy
tCm = 90;

% Total coda window length for Qc and normalization
tWm = 30;

% Parameter for smoothing - must be > 2
nf=8;

% The minimum coda-to-noise energy ratio for the weighted inversion
tresholdnoise = 5;

% Size anomaly for testing: twice(2). Might be set to four(4). The input of
% the checkerboard must be always visually checked visually at the end of
% the process
sizea=2;

% The sped coefficient sets the spectral energy decay of the coda
% wavefield
sped=1.5;

% Number of nodes along x and y
nxc = 10;
nyc = 6;

% How big are markers for stations and events
sz=60;
        
% If the origin time is unknown, you can set a theoretichal velocity for
% the whole area and evaluate it from picking. It must be the velocity of
% the phase you are mapping
vth=3;%in km/s

% Import event origin time and coords of event and station (evst=1)
% from files.
% Ideal is to have both in the SAC header (evst=2) and do it in lat/long.
evst = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TO BE EDITED IF pa>1
% Values of attenuation for testing
hatt =.02;
latt =.001;

% Seconds of time-windows for non-linear inversion and corresponding
% number, set nonlinear=1 to activate;
nonlinear=1;
nW=5;
ntW=tWm/nW;

% Uncertainty on Qc estimation
if nonlinear==0
    % Minimum R-squared for Qc fitting
    RZZ2 = 0.1;
elseif nonlinear==1
    RZZ2 = 2;
    %Grid search
    L1 = 1001;
    m1min=0;
    m1max = 0.01;

    Dmm = (m1max-m1min)/(L1-1);
    m1a = (m1min + Dmm*(0:L1-1))';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name of the file containing the seismograms and rays - in folder "traces"
% creates a list of the waveforms and calculates their number.
% RULE: This is necessary for the option evst = 1, where event and
% station infos are in two separate files. In this case,
% The SAC file-name must have a specific format:
% the first 12 characters are the origin time of the event, as per
% event-location file). Characters 18:20 are the name of station as per
% station file.
%
% For evst = 2, there is no requirement on SAC filenames. In case of pa=3
% the names of rays (if createrays=1) will be identical to file names.

% Get name of files in traces' folder
list = dir(DFolder);
filenames = {list.name}'; %create cell array of file names
filedir = {list.folder}'; %create cell array of file folder
listasac = filenames;
lls = length(listasac);
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

if evst==1
    
   % UTM coordinates of the origin of the model - this must be put at
    % least one cell before the westernmost and southernmost
    % earthquakes/station
    %
    % For Mount Saint Helens
    % the code works with just lat. and long. but the following
    % must be put to zero and the grid must be built accordingly
    originWE=538311;
    originSN=5092338;
    originz=3350;
    
    % Opening the event file. Format is:
    % column (1) = twelve numbers for the origin time of the event
    % (date+time in seconds)
    % column (2) = UTM (WE) or latitude
    % column (3) = UTM (SN) or longitude
    % column (4) = Altitude above sea level in meters
    event=fopen('even.txt');
    name= textscan(event,'%s %f %f %f');
    nameeven=name{1}; even(:,1)=name{2}; even(:,2)=name{3};...
        even(:,3)=name{4};
    fclose(event);

    % Opening the station file. Format is:
    % column (1) = Name of station (3 characters)
    % column (2) = UTM (WE) or latitude
    % column (3) = UTM (SN) or longitude
    % column (4) = Altitude above sea level in meters
    station=fopen('staz.txt');
    name1= textscan(event,'%s %f %f %f');
    namesta=name1{1}; staz(:,1)=name1{2}; staz(:,2)=name1{3};...
        staz(:,3)=name1{4};
    fclose(station);

    numev=length(nameeven); %number of events
    numst=length(namesta); %number of stations

    %Step of the grid for Qc imaging
    stepg=4000;
    degorutm=1;

elseif evst==2
    
    %Latitude and longitude, as usually in Sac Haeder. Stepg in degrees.  
    originWE=20;
    originSN=43;
    originz=0;
    stepg=1;
    degorutm=111;%set to 111 if working with degrees, 1 if in UTM
end

% Build the nodes for peak-delay and Qc imaging -
% they are your parameter-model locations
XY=zeros(floor(nxc)*floor(nyc),2);%pre-define 2D matrix in space

index=0;
for i=1:nxc
    for j=1:nyc
        index=index+1;
        XY(index,1:2)=[originWE+stepg*i originSN+stepg*j];
    end
end
lxy=length(XY(:,1));

% Some figures require moving the point of plotting to half the resolution
halfres=+stepg/2;
x=originWE+halfres:stepg:stepg*(nxc-1)+originWE+halfres;
y=originSN+halfres:stepg:stepg*(nyc-1)+originSN+halfres;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs necessary for the direct-wave CN attenuation tomography
if pa==3
    %1D (1) or 3D (3) velocity model
    dimv=1;
    
    %start time to consider noise - before P-arrival
    sn=5;

    %Rays measured in meters. If the ray length is already in km set um =1;
    um = 1000;

    % This works in [x,y,z], created for UTM WGS84 and origins to zero.
    % In the new version if a 3D velocity model is unavailable a false 3D
    % is created from iasp91
    if dimv==1
        
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
        
    elseif dimv==3
        
        modv=load('modv.txt'); %3D velocity model from text file
        modv(:,5)=0;%Only P-wave info here
        modv(:,1)=modv(:,1)+originWE;
        modv(:,2)=modv(:,2)+originSN;
        
    end
    
    lmod=length(modv(:,3));
    v0= mean(modv(:,4));
    
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
    
    %Creates grid for direct waves and check for zeroes
    %Regular step of the gridD for interpolation - half of step of modv
    
    ixD=floor(passox/resol2x);%numer of x layers,given the step of the gridD
    iyD=floor(passoy/resol2y);%numer of y layers,given the step of the gridD
    izD=floor(passoz/resol2z);%numer of depths layers,given the step of the gridD
    
    gridD=zeros(3,max(passo./resol));
    gridD(1,1:ixD)=originWE:resol2x:originWE+passox-resol2x;
    gridD(2,1:iyD)=originSN:resol2y:originSN+passoy-resol2y;
    gridD(3,1:izD)=-originz:resol2z:-originz+passoz-resol2z;
    
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
    
    %Sets the number of nodes from the origin
    nx=passox/stepg; %how many nodes to change x
    ny=passoy/stepg;

    bs=zeros(nx*ny,2);
    bst=zeros(nx*ny,2);
    lung=zeros(lls,1);
    
    %set createrays to 1 if you want to create files for rays
    %WARNING: This will create A BIG .mat file
    createrays=0;
    
    dtreshold=-3500; %depth treshold for synthetic testing (layers)
end

save('inputs.mat');