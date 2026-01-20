%[text] # MuRAT Multi-Resolution (seismic) Attenuation Tomography
%[text] **SCOPE**: **A code for 3D attenuation-scattering-absorption tomography**
%[text] SYSTEM: The program works on all Mac, Windows and Linux computers where it has been tried.
%[text] MATLAB Version: The code keeps the us of toolboxes to a minimum. However, there are five necessary toolboxes: Signal Processing, Curve Fitting, Image Processing, Mapping and Statistics and Machine Learning Toolboxes. The Parallel Computing Toolbox is highly recommended for speed. Three sample datasets (Mount St. Helens, Romania and Toba) and corresponding inputs are available in this version of MuRAT.
%[text] INSTRUCTIONS: The current release (3.0) works following these steps:
%[text] 1. Download the package at [https://github.com/LucaDeSiena/MuRAT](https://github.com/LucaDeSiena/MuRAT).
%[text] 2. Build your own input file (.mlx) - each field is described in the attached Documentation.pdf and README.txt; however, always use the commented input files as reference. Build your Input files fom those.
%[text] 3. Create your *sac\_data* folder. The utilities Murat\_test.m and Murat\_testAll.m allow you to check that the necessary parameters are in the headers of the SAC files.
%[text] 4. Run this file and select the name of the input file you created beforhand. \
%[text] Author: L. De Siena, December 2025
%[text] ## INPUTS AND CHECKS
fprintf('Checks\n'); %[output:828fc255]

addpath(genpath('Utilities_Matlab'))
addpath('./bin/')

% Ask user for input file
[file, path] = uigetfile({'*.m;*.mlx','MuRAT input files (*.m, *.mlx)'; '*.*','All Files'});
if isequal(file,0)
    error('User selected Cancel');
end
inputFile = fullfile(path, file);
fprintf('User selected %s\n', inputFile); %[output:8e7ea542]

% Run or load input file safely
[~,~,ext] = fileparts(file);
try
    switch lower(ext)
        case '.mlx'
            % Execute live script to set Murat in base or return as struct
            evalc(['run(''' inputFile ''')']); % suppress output
        case '.m'
            % Prefer that the .m returns a struct named Murat; run it
            evalc(['run(''' inputFile ''')']);
        otherwise
            error('Unsupported input file extension: %s', ext);
    end
catch ME
    error('Failed to run/load input file: %s', ME.message);
end
tic;
% Validate and order fields
Murat = Murat_checks(Murat);
Murat.input = orderfields(Murat.input);

% Ensure temp folder exists
tempFolder = fullfile('.', 'temp');
if ~exist(tempFolder,'dir')
    mkdir(tempFolder);
end
save(fullfile(tempFolder,'Murat_checks.mat'), 'Murat');  % default -v7

tCheck = toc;
fprintf('Time spent on checks: %.2f s\n', tCheck);
%%
%[text] ## SEISMIC DATA PROCESSING AND FORWARD MODELLING
fprintf('Data\n');
load(fullfile(tempFolder,'Murat_checks.mat'),'Murat');

% Ask user
useParallel = Murat.input.workers;

% Try to prepare parallel pool only if requested
if ~isempty(useParallel)
    % Check for Parallel Computing Toolbox
    hasPCT = license('test','Distrib_Computing_Toolbox') || license('test','Parallel_Computing_Toolbox');
    if ~hasPCT
        warning('Parallel Computing Toolbox not available. Falling back to sequential.');
        useParallel = false;
    else
        % Start pool if none exists
        pool = gcp('nocreate');
        if isempty(pool)
            try
                pool = useParallel; % optionally specify NumWorkers
            catch ME
                warning('Failed to start parpool: %s. Falling back to sequential.');
                useParallel = false;
            end
        end
    end
end

tic;

if ~isempty(useParallel)
    Murat = Murat_dataParallelized(Murat);
else
    Murat = Murat_data(Murat);
end
save(fullfile(tempFolder,'Murat_forward.mat'), 'Murat');
tData = toc;
fprintf('Time spent on data processing: %.2f s\n', tData);
%%
%[text] ## TOMOGRAPHIC INVERSIONS
tic;
load(fullfile(tempFolder,'Murat_forward.mat'),'Murat');
fprintf('Inversion\n');
Murat = Murat_inversion(Murat);
save(fullfile(tempFolder,'Murat_inverse.mat'), 'Murat');
tInv = toc;
fprintf('Time spent on inversion: %.2f s\n', tInv);
%%
%[text] ## CREATING PLOTS and CLEANING UP
tic;
load(fullfile(tempFolder,'Murat_inverse.mat'),'Murat');
fprintf('Plot\n');
Murat = Murat_plot(Murat);

% Ensure output folder exists and save final Murat (use -v7.3 if large)
outFolder = Murat.input.label;
if ~exist(outFolder,'dir')
    mkdir(outFolder);
end
save(fullfile(outFolder,'Murat.mat'),'Murat');  % change to -v7.3 if >2GB

tPlot = toc;
fprintf('Time spent on plotting: %.2f s\n', tPlot);

% Clean up only temporary large variables
clearvars -except Murat tCheck tData tInv tPlot

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:828fc255]
%   data: {"dataType":"text","outputData":{"text":"Checks\n","truncated":false}}
%---
%[output:8e7ea542]
%   data: {"dataType":"text","outputData":{"text":"User selected \/Users\/lucadesiena\/Documents\/Documents - Luca’s MacBook Pro (2)\/MATLAB\/MuRAT_December_2025\/Murat_inputToba.m\n","truncated":false}}
%---
