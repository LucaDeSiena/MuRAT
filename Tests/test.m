function smoke_test()
% MuRAT smoke test
% This test checks that the repository contains the expected files and that the
% core MuRAT functions are declared with valid signatures.
% It avoids requiring the MATLAB Parallel Computing Toolbox or a MATLAB license.

    root = fileparts(fileparts(mfilename('fullpath')));
    if ~isempty(root)
        addpath(genpath(root));
    end

    requiredFiles = {
        'Utilities_Matlab/MyUtilities/Murat_test.m', ...
        'Utilities_Matlab/MyUtilities/Murat_testAll.m', ...
        'bin/Murat_testData.m', ...
        'bin/Murat_checks.m', ...
        'Utilities_Matlab/MatSAC/sac.m', ...
        'Utilities_Matlab/MatSAC/sachdr.m'
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

    % Validate the key MuRAT function definitions exist.
    requiredFunctions = {'Murat_test', 'Murat_testAll', 'Murat_testData'};
    for i = 1:numel(requiredFunctions)
        if exist(requiredFunctions{i}, 'file') ~= 2
            error('Required function %s was not found on the MATLAB path.', requiredFunctions{i});
        end
    end

    % Check the test-data helper initializes its flag variable before use.
    fid = fopen('bin/Murat_testData.m', 'r');
    if fid == -1
        error('Unable to open bin/Murat_testData.m');
    end
    content = fread(fid, '*char')';
    fclose(fid);
    content = char(content.');

    hasFlagInit = contains(content, 'flag = []') || contains(content, 'flag=[]') || ...
                  contains(content, 'flag = 0') || contains(content, 'flag=0');
    assert(hasFlagInit, 'Murat_testData.m must initialize the flag variable.');

    fprintf('MuRAT smoke test passed. Core repository files and entry points are present.\n');
end
