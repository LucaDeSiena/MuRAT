function test()
    % MuRAT smoke test for CI.
    % This test checks the repository structure and core MuRAT entries
    % without needing the Parallel Computing Toolbox or a MATLAB license.

    root = fileparts(fileparts(mfilename('fullpath')));
    if isempty(root)
        root = pwd;
    end

    addpath(root);
    addpath(fullfile(root, 'bin'));
    addpath(fullfile(root, 'Utilities_Matlab'));
    addpath(fullfile(root, 'Utilities_Matlab', 'MatSAC'));
    addpath(fullfile(root, 'Utilities_Matlab', 'MyUtilities'));

    requiredFiles = {
        fullfile(root, 'Utilities_Matlab', 'MyUtilities', 'Murat_test.m');
        fullfile(root, 'Utilities_Matlab', 'MyUtilities', 'Murat_testAll.m');
        fullfile(root, 'bin', 'Murat_testData.m');
        fullfile(root, 'bin', 'Murat_checks.m');
        fullfile(root, 'Utilities_Matlab', 'MatSAC', 'sac.m');
        fullfile(root, 'Utilities_Matlab', 'MatSAC', 'sachdr.m')
    };

    missingFiles = {};
    for i = 1:numel(requiredFiles)
        fname = requiredFiles{i};
        if exist(fname, 'file') ~= 2
            missingFiles{end+1} = fname;
            fprintf('  X MISSING: %s\n', fname);
        else
            fprintf('  + Found: %s\n', fname);
        end
    end

    if ~isempty(missingFiles)
        error('Missing required files');
    end

    requiredFunctions = {'Murat_test', 'Murat_testAll', 'Murat_testData'};
    for i = 1:numel(requiredFunctions)
        fname = requiredFunctions{i};
        if exist(fname, 'file') ~= 2
            error('Required function not found: %s', fname);
        end
    end

    fid = fopen(fullfile(root, 'bin', 'Murat_testData.m'), 'r');
    if fid == -1
        error('Unable to open bin/Murat_testData.m');
    end
    content = fread(fid, '*char');
    fclose(fid);
    content = char(content.');

    hasFlagReturn = ~isempty(strfind(content, 'function [muratHeader,flag]'));
    hasFlagInit = ~isempty(strfind(content, 'flag = []')) || ...
                 ~isempty(strfind(content, 'flag=[]')) || ...
                 ~isempty(strfind(content, 'flag = 0')) || ...
                 ~isempty(strfind(content, 'flag=0'));

    if ~hasFlagReturn
        error('Murat_testData.m must return [muratHeader, flag]');
    end
    if ~hasFlagInit
        error('Murat_testData.m must initialize flag variable');
    end

    fprintf('MuRAT smoke test passed.\n');
end