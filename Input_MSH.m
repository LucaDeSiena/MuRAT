%% INPUTS
addpath('/Users/lucadesiena/Documents/Utilities_Matlab')

clear
close all
clc

disp('Input Section')

% The following prompts are included TO BE  EDITED before running!

%Which analysis do you want to perform?
%Pick delay and Qc without kernels: pa=1
%Pick delay and Qc with kernels: pa=2
%Pick delay, kernel-Qc and P/S wave attenuation with the CN method: pa=3
pa=3;

%Folder containing data
DFolder = '/Users/lucadesiena/Documents/Datasets/MSH_2000_2003/*.SAC';

%Creates folders and paths to store results and figures
FLabel = 'MSH';
FPath = '/Users/lucadesiena/Documents/MATLAB/MURAT5_3_6';

%mkdir(FLabel)

% Storing inversion parameters
inputinv = 'inversion.mat';
inputdata = 'data.mat';

%Output figure format 
fformat='jpeg';

%Figure visibility - pick 'on' or 'off'
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
cf = 6;

% Maximum window to pick pick-delays in seconds
maxtpde = 10;

% Minimum peak delay considering scattering in the area and frequency
mintpde = 1/cf;

% Lapse time for the start of the window used to measure and calculate the
% normalization energy
tCm = 30;

% Total coda window length for Qc and normalization
tWm = 15;

% Parameter for smoothing - must be > 2
nf=8;

% The minimum coda-to-noise energy ratio for the weighted inversion
tresholdnoise = 10;

% Size anomaly for testing: twice(2). Might be set to four(4). The input of
% the checkerboard must be always visually checked visually at the end of
% the process
sizea=2;

% The sped coefficient sets the spectral energy decay of the coda
% wavefield
sped=0.5;

% Number of nodes along x and y - must be even for checkerboard to work!
nxc = 10;
nyc = 12;

% Size of markers for stations and events on maps
sz=60;
        
% If the origin time is unknown, you can set a theoretichal velocity for
% the whole area and evaluate it from picking. It must be the velocity of
% the phase you are mapping. 
vth=3;%in km/s

% Import event origin time and coords of event and station (evst=1)
% from files.
% Ideal is to have both in the SAC header (evst=2) and do it in lat/long.
evst = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TO BE EDITED IF pa>1
% Values of attenuation for testing in checkerboard
hatt =.02;
latt =.001;

% Seconds of time-windows for non-linear inversion and corresponding
% number, set nonlinear=1 to activate;
nonlinear=1;
nW=3;
ntW=tWm/nW;

% Uncertainty on Qc estimation
if nonlinear==0
    % Minimum R-squared for Qc fitting
    RZZ2 = 0.1;
elseif nonlinear==1
    RZZ2 = 5;
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

elseif evst==2
    
    %Latitude and longitude, as usually in Sac Haeder. Stepg in degrees.
    
    originWE=20;
    originSN=43;
    originz=0;
    um=1;
    stepg=1;% degrees
    
end

% Build the nodes for peak-delay and Qc imaging -
% they are your parameter-model locations
XY=zeros(nxc*nyc,2);%pre-define 2D matrix in space

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
    
    %start time to consider noise - before P-arrival
    sn=5;

    %Rays measured in meters. If the ray length is already in km set um =1;
    um = 1000;

    %This works in [x,y,z] or better UTM WGS84, setting origins to zero
    modv=load('modv.txt');
    modv(:,5)=0;%Only P-wave info here
    modv(:,1)=modv(:,1)+originWE; 
    modv(:,2)=modv(:,2)+originSN; 
    
    lmod=length(modv(:,3));
    v0= mean(modv(:,4));
    resol2 = abs(modv(2,3)-modv(1,3))/2;
    
    %Steps of the velocity models
    passox=max(modv(:,1))-min(modv(:,1));
    passoy=max(modv(:,2))-min(modv(:,2));
    passoz=max(modv(:,3))-min(modv(:,3));
    
    passo=[passox passoy passoz];
    
    %Creates grid for direct waves and check for zeroes
    gridD=zeros(3,max(passo/resol2));
    gridD(1,1:passox/resol2)=originWE:resol2:originWE+passox-resol2;
    gridD(2,1:passoy/resol2)=originSN:resol2:originSN+passoy-resol2;
    gridD(3,1:passoz/resol2)=-originz:resol2:-originz+passoz-resol2;
    
    %Regular step of the gridD for interpolation - half of step of modv
    
    ixD=passox/resol2;%numer of x layers,given the step of the gridD
    iyD=passoy/resol2;%numer of y layers,given the step of the gridD
    izD=passoz/resol2;%numer of depths layers,given the step of the gridD
    
    % gridD dimensions
    pvel = zeros(ixD,iyD,izD);
    
    %NUMBER OF X, Y, AND Z LAYERS
    for k=1:izD
        m1=[modv(:,1:3) modv(:,PorS+2)];
        index=0;
        for i=1:ixD
            for j=1:iyD
                index=index+1;
                pvel(i,j,k) = m1(index,4);
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
    
    WEi=565000;
    SNi=5115000;
    zi=-2000;
end

save('inputs.mat');