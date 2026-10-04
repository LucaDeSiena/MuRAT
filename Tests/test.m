function test()
% MuRAT smoke test for CI.
% This test checks that the repository contains the expected MuRAT files and
% that the key function entry points are discoverable on the MATLAB path.

    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    if isempty(repoRoot)
        repoRoot = pwd;
    end

    addpath(repoRoot);
    addpath(fullfile(repoRoot, 'bin'));
    addpath(fullfile(repoRoot, 'Utilities_Matlab'));
    addpath(fullfile(repoRoot, 'Utilities_Matlab', 'MatSAC'));
    addpath(fullfile(repoRoot, 'Utilities_Matlab', 'MyUtilities'));

    requiredFiles = {
        fullfile(repoRoot, 'Utilities_Matlab', 'MyUtilities', 'Murat_test.m'), ...
        fullfile(repoRoot, 'Utilities_Matlab', 'MyUtilities', 'Murat_testAll.m'), ...
        fullfile(repoRoot, 'bin', 'Murat_testData.m'), ...
        fullfile(repoRoot, 'bin', 'Murat_checks.m'), ...
        fullfile(repoRoot, 'Utilities_Matlab', 'MatSAC', 'sac.m'), ...
        fullfile(repoRoot, 'Utilities_Matlab', 'MatSAC', 'sachdr.m')
    };

    missing = {};
    for i = 1:numel(requiredFiles)
        if exist(requiredFiles{i}, 'file') ~= 2
            missing{end+1} = requiredFiles{i};
        end
    end

    if ~isempty(missing)
        error('Missing required MuRAT files: %s', strjoin(missing, ', '));
    end

    requiredFunctions = {'Murat_test', 'Murat_testAll', 'Murat_testData'};
    for i = 1:numel(requiredFunctions)
        if exist(requiredFunctions{i}, 'file') ~= 2
            error('Required function %s was not found on the MATLAB path.', requiredFunctions{i});
        end
    end

    fid = fopen(fullfile(repoRoot, 'bin', 'Murat_testData.m'), 'r');
    if fid == -1
        error('Unable to open bin/Murat_testData.m');
    end
    content = fread(fid, '*char');
    fclose(fid);
    content = char(content.');

    hasFlagInit = contains(content, 'flag = []') || contains(content, 'flag=[]') || ...
                  contains(content, 'flag = 0') || contains(content, 'flag=0');
    assert(hasFlagInit, 'Murat_testData.m must initialize the flag variable.');

    fprintf('MuRAT smoke test passed. Core repository files and entry points are present.\n');
end