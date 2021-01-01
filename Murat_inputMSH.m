%% Hard coded data in .m spreadsheet - 3D Analysis
% This is an input file for the program Multi-Resolution Attenuation
% Tomography (MuRAT). It refers to the following area:
%
%   MOUNT ST HELENS VOLCANO
%
%   EVERYTHING MARKED BY 'E' IS TO/CAN BE EDITED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PATHS AND FIGURES - GENERAL FIELD AND FIGURES
% INPUTS in this section set the directory space where you want to work as
% well as the characteristics of the figure files you output.

% Working directory
Murat.input.workingDirectory        =	'./';%'E'

% Folder containing data
Murat.input.dataDirectory           =   './sac_MSH/*.sac';%'E'

% Name of folder to store results and figures
Murat.input.label                   =   'MSH';%'E'

% Choose if figures related to data processing are to be created.
% If set to 1 it will output the test for stability of peak delay and Qc
Murat.input.data                    =   1;

% Choose if figures related to the velocity model are to be created.
% If set to 1 it will output the existing velocity model, rays, and Qc
% sensitivity kernels 
Murat.input.geometry                =	1;

% Choose if figures related to the the inversion results are to be created.
% If set to 1 it will output the l-curves and will allow to input the
% damping, otherwise it must be set beforehand.
Murat.input.inversion               =	1;%E

% Choose if figures related to the checkerboard test are to be created.
% If set to 1 it will output the input and output of the checkerboards for 
% Q and Qc.
Murat.input.checkerboard            =   1;

% Choose if figures related to the spike test are to be created.
% If set to 1 it will output the input and output of the spike tests for 
% Q and Qc.
Murat.input.spike                   =	1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DATA DRIVEN CHOICES - WAVEFORMS
% INPUTS in this section set all the variables required by data processing.
% This includes data choises, as setting the name of the variables in SAC
% and all the attributes that are needed for the three kinds of imaging.


% ----
% WAVEFORMS
% Name of the SAC variables where zero time and P/S pickings are saved.
Murat.input.originTime              =   [];%'E'
Murat.input.PTime                   =   'SAChdr.times.a';%'E'
Murat.input.STime                   =	[];%'E'

% Choose what phase you are interested in, P- (2) or S- (3).
Murat.input.POrS                    =	2; %'E'

% Central frequency (Hz) - set it according to your spectrograms.
Murat.input.centralFrequency        =	6;%'E'

% Work with 1 vertical (1) or 2 horizontal (2) recordings, or with the
% three components(3). Order MUST BE: WE, SN, Z!
Murat.input.components              =	1;%'E'

% Parameter for smoothing - must be > 2 and an integer.
Murat.input.smoothing               =	8;%'E'

% ----
% PEAK DELAY AND DIRECT WAVE
% Minimum peak delay considering scattering in the area and frequency
Murat.input.minimumPeakDelay        =	0.1;%'E'

% Maximum window to pick pick-delays in seconds
Murat.input.maximumPeakDelay        =	7;%'E'

% DEVELOPMENT: Order of the modified Bessel function of second kind,
% describing the von-Kármán type random medium
% Murat.input.vonKarman              =   0.5;
% 

% The spectral energy decay of the coda wavefield, either surface (0.5) or
% body waves (1)
Murat.input.spectralDecay           =	0.5;%'E'

% Length of the window used to measure noise and direct wave intensity.
Murat.input.bodyWindow              =	1;

% Start of the window used to measure noise from origin of record.
Murat.input.startNoise              =	5;

% The minimum coda-to-noise energy ratio for the weighted inversion.
Murat.input.tresholdNoise           =	3;

% ----
% Qc
% Lapse time for the start of the window used to measure and calculate the
% normalization energy
Murat.input.startLapseTime          =   15;%'E'

% Total coda window length for Qc and normalization
Murat.input.codaWindow              =   15;%'E'

% Best fitting Qc to model intensity data measured on smaller windows
% with a grid search. Set Murat.input.nonLinear == 1 to activate.
% Classical mode is ==0.
Murat.input.nonLinear               =   0;%'E'

if Murat.input.nonLinear == 0
    
    % Minimum R-squared for Qc fitting
    Murat.input.fitTreshold         =	0.1;%'E'
    
elseif Murat.input.nonLinear == 1
    
    % How many times you divide the total coda window 
    Murat.input.fitNumberWindows    =   3;%'E'
    
    % Length of smaller time windows to compute coda intensity    
    Murat.input.fitLengthWindows    =   5;%'E'
    
    % minimum inverse Qc allowed
    Murat.input.minimumQcGrid       =   0;%E
    
    % maximum inverse Qc allowed
    Murat.input.maximumQcGrid       =   0.01;%E
    
    % total number of Qc in the interval
    Murat.input.totalQcNumber       =   1001;%E
    
end

% Treshold to reduce computational time for Pacheco-Snieder kernels.
% It divides the inversion grid by the treshold.
% Set it if it takes too long to run the data computation.
Murat.input.kernelTreshold          =	2;%'E'

%----
%Geometry and velocity
% Coordinates of the origin of the model in lat/lon  - must be put
% (1) at least one cell before the westernmost and southernmost
% earthquakes/station for positive latitudes and longitudes
%
% (2) at least one cell after the easternmost and northermost
% earthquakes/station for positive latitudes and longitudes
% if the origin is given in UTM then convert to degrees: use utm2deg
% Murat.input.origin                 =   [538500 5092000 3350];%'E'
% Murat.input.end                    =   [586500 5138000 -22650];%'E'

% Give latitude, longitude and depth (m) of the origin and end of the model
Murat.input.origin                  =   [45.9805 -122.5030 3350];%'E'
Murat.input.end                     =   [46.3900 -122.0000 -22650];%'E'

% Number of nodes in the three directions.
Murat.input.gridX                   =   10;%'E'
Murat.input.gridY                   =   12;%'E'
Murat.input.gridZ                   =   6;%'E'

%Choose if you want to use none (0) or a 3D (1) velocity model
Murat.input.availableVelocity       =	1;%'E'

% Set a theoretichal velocity for % the whole area if you have no
% available velocity model or have no info of origin time.
% It must be the velocity of the phase you are mapping in km/s.
Murat.input.averageVelocityP        =   6;%'E'
Murat.input.averageVelocityS        =   3;%'E'

%name of the velocity model - if option (1)
if Murat.input.availableVelocity == 1
    Murat.input.namev               =   'modv.txt';%'E'
end

% Only needed if you want to import external files with locations of
% events and stations - in this case set to 1.
Murat.input.importLocation          =   1;%DEPRECATED but still available.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RESULTS - INCLUDES INVERSION, TESTING AND PLOTTING 
if Murat.input.inversion == 0
    Murat.input.lcurveQc            =   0.15;%E
    Murat.input.lcurveQ             =   0.03;%E
end

% Testing - checkerboard
% Size anomaly for testing: twice(2) or four times node spacing.
% The input of the checkerboard must be  checked visually at the end of
% the process    
Murat.input.sizeCheck               =   2;%'E'

% Values of attenuation for testing
Murat.input.highCheck               =   0.02;%'E'
Murat.input.lowCheck                =   0.001;%'E'

% Testing - checkerboard
% Set the origin ad end locations for the spike test
Murat.input.spikeLocationOrigin     =   [46.05 -122.40  0];%'E'
Murat.input.spikeLocationEnd        =   [46.12 -122.30 -6000];%'E'

%Set the value of the spike compared to the average
Murat.input.spikeValue              =   0.02;%'E'

% Figure format - 'jpeg' (fast) or 'tiff' (for publication)
Murat.input.format                  =   'jpeg';%'E'

% Figure visibility during computation - pick 'on' or 'off'
Murat.input.visibility              =   'on';%'E'

% Sections to slice plots
Murat.input.sections                =	[46.12 -122.20 -1000];

% How big are markers for stations and events
Murat.input.sizeMarker              =   60;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
