%[text] # MuRAT Multi-Resolution (seismic) Attenuation Tomography
%[text] **SCOPE**: **A code for 3D attenuation-scattering-absorption tomography**
%[text] SYSTEM: The program works on all Mac, Windows and Linux computers where it has been tried.
%[text] MATLAB Version: Matlab 2025b
%[text] TOOLBOXES:
%[text] - Curve Fitting Toolbox
%[text] - Mapping Toolbox
%[text] - Optimization and Global Optimization Toolboxes
%[text] - Signal Processing Toolbox.
%[text] - Parallel Computing Toolbox is highly recommended for speed. \
%[text] TEMPLATES: Three sample datasets (Mount St. Helens, Romania and Toba) and corresponding inputs are available in the root MuRAT folder.
%[text] INSTRUCTIONS: The current release (4.0) works following these steps:
%[text] 1. Download or branch the package at [https://github.com/LucaDeSiena/MuRAT](https://github.com/LucaDeSiena/MuRAT).
%[text] 2. Build your own input file (.m or .mlx) - each field is described in the attached Documentation.pdf and README.txt.
%[text] 3. Create your *sac\_data* folder. The codes in the *Utilities\_Matlab* folder, especially those in the *PrePostProcessing* subfolder, can be useful to check data and outputs, especially the necessary parameters are in the headers of the SAC files.
%[text] 4. Run this file and select the name of the input file you created beforhand. \
%[text] Author: L. De Siena, January 2026
%[text] ## INPUTS AND CHECKS
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

% Ask user for input file
[file, path]= uigetfile({'*.m;*.mlx','MuRAT input files (*.m, *.mlx)'; '*.*','All Files'});
if isequal(file,0)
    error('No input file selected.');
end
inputFile   = fullfile(path, file);
fprintf('Using input file %s\n', inputFile);

% Run or load input file safely
[~,~,ext]   = fileparts(file);
switch lower(ext)
    case '.mlx'
        run(inputFile);
    case '.m'
        run(inputFile);
    otherwise
        error('Unsupported input file extension: %s', ext);
end
tStart      = tic;

% --- Validation and ordering ---
Murat       =   Murat_checks(Murat);
Murat.input =   orderfields(Murat.input);

% --- Ensure temp folder exists ---
tempFolder  = fullfile(pwd,'temp');
if ~exist(tempFolder,'dir'), mkdir(tempFolder); end
save(fullfile(tempFolder,'Murat_checks.mat'), 'Murat','-v7');

tCheck = toc(tStart);
fprintf('Time spent on checks: %.2f s\n', tCheck);
%%
%[text] ## SEISMIC DATA PROCESSING AND FORWARD MODELLING
load(fullfile(tempFolder,'Murat_checks.mat'),'Murat');

% Ask user
useParallel     =   Murat.input.workers;

% Try to prepare parallel pool only if requested
if ~isempty(useParallel)
    % Check for Parallel Computing Toolbox
    hasPCT      =   license('test','Distrib_Computing_Toolbox') || license('test','Parallel_Computing_Toolbox');
    if ~hasPCT
        warning('Parallel Computing Toolbox not available. Falling back to sequential.');
        useParallel =   false;
    else
        % Start pool if none exists
        pool        =   gcp('nocreate');
        if isempty(pool)
            try
                pool= useParallel; % optionally specify NumWorkers
            catch ME
                warning('Failed to start parpool: %s. Falling back to sequential.');
                useParallel = false;
            end
        end
    end
end

tDataStart = tic;

if useParallel
    Murat   =   Murat_dataParallelized(Murat);
else
    Murat   =   Murat_data(Murat);
end
save(fullfile(tempFolder,'Murat_forward.mat'), 'Murat','-v7');
tData = toc(tDataStart);
fprintf('Time spent on data processing: %.2f s\n', tData);
%%
%[text] ## TOMOGRAPHIC INVERSIONS
tInvStart   =   tic;
load(fullfile(tempFolder,'Murat_forward.mat'),'Murat');
Murat       =   Murat_inversion(Murat);
% Save with -v7.3 if large
info        =   whos('Murat');
if info.bytes > 2e9
    save(fullfile(tempFolder,'Murat_inverse.mat'),'Murat','-v7.3');
else
    save(fullfile(tempFolder,'Murat_inverse.mat'),'Murat','-v7');
end
tInv        =   toc(tInvStart);
fprintf('Inversion completed (%.2f s)\n', tInv);
%%
%[text] ## PLOTTING and CLEANING UP
tPlotStart  =   tic;
Murat       =   Murat_plot(Murat);

outFolder   =   Murat.input.label;

info        =   whos('Murat'); % re-check size
if info.bytes > 2e9
    save(fullfile(outFolder,'Murat.mat'),'Murat','-v7.3');
else
    save(fullfile(outFolder,'Murat.mat'),'Murat','-v7');
end

tPlot       =   toc(tPlotStart);
fprintf('Plotting completed (%.2f s)\n', tPlot);

totalTime   =   toc(tStart);
fprintf('Total run time: %.2f s\n', totalTime);

clearvars -except Murat tCheck tData tInv tPlot totalTime

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
