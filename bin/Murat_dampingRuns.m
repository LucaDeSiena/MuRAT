function [Murat, dampingRuns, tInv, tPlot] = Murat_dampingRuns(Murat, tempFolder, phase)
% MURAT_DAMPINGRUNS Invert and plot one model or a list of damping runs.
%
% Call with phase = 'prepare' before Murat_checks to select the first
% damping folder for forward-stage output; omit phase for inversion/plotting.
% Inputs:
%   dampingSweepValues: [] for normal processing; otherwise a vector of
%       nonnegative damping values, one separate run per value. Each value
%       applies to BOTH Qc and Q at all existing central frequencies.
% New sweeps archive existing matching folders as _old, _old_2, etc.
%
% With an empty/missing sweep, dampingValueQc and dampingValueQ retain their
% normal per-frequency meanings. A sweep overrides those two fields only.
% Smoothing, frequencies, and all other analysis settings are preserved.
%
% Outputs: Murat is the LAST completed model, dampingRuns is a lightweight
% summary, and tInv/tPlot are cumulative inversion/plotting times in seconds.
% Every completed run saves its full Murat structure separately on disk.

if nargin < 3, phase = 'run'; end
if ~ismember(phase, {'prepare','run'})
    error('MuRAT:dampingPhase', 'Phase must be prepare or run.');
end
dampingRuns = []; tInv = 0; tPlot = 0;
if nargin < 2
    tempFolder = fullfile(pwd, 'temp');
end
baseMurat = Murat;
sourceLabel = char(Murat.input.label);
baseLabel = sourceLabel;
if isfield(Murat.input, 'dampingBaseLabel')
    baseLabel = char(Murat.input.dampingBaseLabel);
end
frequencies = Murat.input.centralFrequency;
values = [];
if isfield(Murat.input, 'dampingSweepValues')
    values = Murat.input.dampingSweepValues;
end
isSweep = ~isempty(values);
if ~isSweep, baseLabel = sourceLabel; end

if isSweep
    validateattributes(values, {'numeric'}, ...
        {'vector','real','finite','nonnegative'}, mfilename, 'dampingSweepValues');
    methods = {'Tikhonov', 'Particle', 'Annealing', 'Global'};
    match = find(strcmpi(strtrim(char(Murat.input.inversionMethod)), methods), 1);
    if isempty(match)
        error('MuRAT:dampingSweepMethod', ...
            'Supported inversion methods: Tikhonov, Particle, Annealing, Global.');
    end
    Murat.input.inversionMethod = methods{match};
    baseMurat.input.inversionMethod = methods{match};
    values = double(values(:).');
    if numel(unique(values)) ~= numel(values)
        error('MuRAT:duplicateDamping', 'Use each damping value only once.');
    end
    nRuns = numel(values);
else
    nRuns = 1;
end

names = cell(1, nRuns);
folders = cell(1, nRuns);
for i = 1:nRuns
    if isSweep
        % Preserve familiar names such as damp_0p005, with enough precision
        % to distinguish closely spaced double-precision input values.
        names{i} = ['damp_' strrep(strrep(sprintf('%.15g', values(i)), ...
            '.', 'p'), '-', 'm')];
        folders{i} = [baseLabel '_' names{i}];
    else
        names{i} = 'single';
        folders{i} = baseLabel;
    end
end
if numel(unique(folders)) ~= nRuns
    error('MuRAT:dampingNameCollision', ...
        'Damping values produce identical folder names; choose distinct values.');
end
% Check ALL destinations before starting any inversion.
for i = 1:nRuns
    if isfile(folders{i})
        error('MuRAT:outputIsFile', 'Output path is a file: %s', folders{i});
    end

end

% Select the first requested damping folder BEFORE checks create outputs.
if strcmp(phase, 'prepare')
    if isSweep
        % Archive all matching runs before checks create the first folder.
        for i = 1:nRuns
            if isfolder(folders{i})
                archive = [folders{i} '_old'];
                archiveIndex = 2;
                while isfolder(archive) || isfile(archive)
                    archive = sprintf('%s_old_%d', folders{i}, archiveIndex);
                    archiveIndex = archiveIndex + 1;
                end
                [ok, message] = movefile(folders{i}, archive);
                if ~ok
                    error('MuRAT:dampingArchiveFailed', ...
                        'Could not archive %s: %s', folders{i}, message);
                end
                fprintf('Archived previous run: %s -> %s\n', folders{i}, archive);
            end
        end
        Murat.input.dampingBaseLabel = baseLabel;
        Murat.input.label = folders{1};
    end
    return
end

% Include every directory used by checks, inversions, and plotting.
subdirs = {'', 'TXT', 'MAT', 'FIGURES', 'matfiles', 'Rays', 'Kernels', ...
    'Tests', 'Tests/PeakDelay', 'Tests/Qc', 'Tests/Q', 'Tests/LCurve', ...
    'Tests/Clustering', 'Results', 'Results/PeakDelay', 'Results/Qc', ...
    'Results/Q', 'Results/Parameter', 'Results/Residuals', 'Results/L-curve', ...
    'Checkerboard', 'Checkerboard/PD', 'Checkerboard/Qc', 'Checkerboard/Q', ...
    'Spike', 'Spike/Qc', 'Spike/Q'};
dampingRuns = repmat(struct('runName', '', 'label', '', 'dampingValueQc', [], ...
    'dampingValueQ', [], 'cacheFile', '', 'inversionSeconds', 0, ...
    'plottingSeconds', 0), nRuns, 1);
tInv = 0;
tPlot = 0;

for i = 1:nRuns
    % Reset to the SAME observed data, kernels, and forward model each time.
    Murat = baseMurat;
    Murat.input.label = folders{i};
    if isSweep
        Murat.input.dampingRunName = names{i};
        Murat.input.dampingValueQc = values(i) * ones(size(frequencies));
        Murat.input.dampingValueQ = values(i) * ones(size(frequencies));
    end
    for j = 1:numel(subdirs)
        destination = fullfile(folders{i}, subdirs{j});
        if ~isfolder(destination)
            mkdir(destination);
        end
    end
    if isSweep && ~strcmp(sourceLabel, folders{i})
        % Copy only forward-stage diagnostics, not previous inverse models.
        for diagnostic = {'Kernels', 'Rays', 'Tests'}
            source = fullfile(sourceLabel, diagnostic{1});
            target = fullfile(folders{i}, diagnostic{1});
            if isfolder(source)
                items = dir(source);
                items = items(~ismember({items.name}, {'.','..'}));
                for j = 1:numel(items)
                    copyfile(fullfile(source, items(j).name), ...
                        fullfile(target, items(j).name));
                end
            end
        end
    end

    fprintf('\nMuRAT run %d/%d: %s\n', i, nRuns, folders{i});
    timer = tic;
    Murat = Murat_inversion(Murat);
    invSeconds = toc(timer);
    tInv = tInv + invSeconds;
    fprintf('Inversion completed (%.2f s)\n', invSeconds);
    if ~isSweep
        if ~isfolder(tempFolder), mkdir(tempFolder); end
        saveMurat(fullfile(tempFolder, 'Murat_inverse.mat'), Murat);
    end

    timer = tic;
    Murat = Murat_plot(Murat);
    plotSeconds = toc(timer);
    tPlot = tPlot + plotSeconds;
    fprintf('Plotting completed (%.2f s)\n', plotSeconds);
    if isSweep
        cacheFile = fullfile(folders{i}, ['Murat_' names{i} '.mat']);
    else
        cacheFile = fullfile(folders{i}, 'Murat.mat');
    end
    saveMurat(cacheFile, Murat);
    dampingRuns(i) = struct('runName', names{i}, 'label', folders{i}, ...
        'dampingValueQc', Murat.input.dampingValueQc, ...
        'dampingValueQ', Murat.input.dampingValueQ, 'cacheFile', cacheFile, ...
        'inversionSeconds', invSeconds, 'plottingSeconds', plotSeconds);
    fprintf('Saved completed run: %s\n', cacheFile);
end
end

function saveMurat(filename, Murat)
% Use the same size-based MAT format selection as the original main script.
info = whos('Murat');
if info.bytes > 2e9
    save(filename, 'Murat', '-v7.3');
else
    save(filename, 'Murat', '-v7');
end
end
