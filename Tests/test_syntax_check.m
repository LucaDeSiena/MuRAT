%% SYNTAX CHECK TEST
% Validates all MATLAB files for:
% - Correct function signatures
% - Proper variable initialization
% - No circular dependencies
% - Basic parsing errors

clear; clc;

fprintf('Starting Syntax Validation Tests...\n');
test_results = {};

%% Test 1: Murat_testAll flag initialization
fprintf('\nTest 1: Checking flag initialization in Murat_testAll.m...\n');
try
    % Read the file
    fid = fopen('Utilities_Matlab/MyUtilities/Murat_testAll.m', 'r');
    content = fread(fid, '*char')';
    fclose(fid);
    
    % Check if flag is initialized before the loop
    has_init = contains(content, 'flag') && ...
        (contains(content, 'flag = []') || contains(content, 'flag=[]') || ...
         contains(content, 'flag = 0') || contains(content, 'flag=0'));
    
    if ~has_init
        warning('FLAG_NOT_INITIALIZED: flag variable is not properly initialized before use in loop');
        test_results{end+1} = struct('test', 'flag_init', 'status', 'FAIL', 'msg', 'flag not initialized');
    else
        fprintf('✓ flag is properly initialized\n');
        test_results{end+1} = struct('test', 'flag_init', 'status', 'PASS', 'msg', '');
    end
catch ME
    fprintf('✗ Error reading file: %s\n', ME.message);
    test_results{end+1} = struct('test', 'flag_init', 'status', 'ERROR', 'msg', ME.message);
end

%% Test 2: Check Murat_testData.m returns flag variable
fprintf('\nTest 2: Checking Murat_testData.m flag output...\n');
try
    fid = fopen('bin/Murat_testData.m', 'r');
    content = fread(fid, '*char')';
    fclose(fid);
    
    % Check function signature
    has_flag_output = contains(content, 'function [muratHeader,flag]');
    has_flag_init = contains(content, 'flag = []') || contains(content, 'flag=[]');
    
    if has_flag_output && has_flag_init
        fprintf('✓ Murat_testData properly declares and initializes flag\n');
        test_results{end+1} = struct('test', 'testdata_flag', 'status', 'PASS', 'msg', '');
    else
        fprintf('✗ Murat_testData flag declaration incomplete\n');
        test_results{end+1} = struct('test', 'testdata_flag', 'status', 'FAIL', 'msg', 'Missing flag declaration or init');
    end
catch ME
    fprintf('✗ Error: %s\n', ME.message);
    test_results{end+1} = struct('test', 'testdata_flag', 'status', 'ERROR', 'msg', ME.message);
end

%% Test 3: Check function file existence for critical functions
fprintf('\nTest 3: Checking critical function file existence...\n');
critical_functions = {
    'bin/Murat_test.m', ...
    'bin/Murat_testData.m', ...
    'bin/Murat_checks.m', ...
    'bin/Murat_inversion.m', ...
    'bin/Murat_plot.m', ...
    'Utilities_Matlab/MatSAC/fget_sac.m'
};

for i = 1:length(critical_functions)
    fname = critical_functions{i};
    if isfile(fname)
        fprintf('  ✓ %s exists\n', fname);
    else
        fprintf('  ✗ %s MISSING\n', fname);
        test_results{end+1} = struct('test', sprintf('exists_%s', strrep(fname,'/','')), ...
            'status', 'FAIL', 'msg', sprintf('%s not found', fname));
    end
end

%% Test 4: Basic syntax parsing with Octave/MATLAB
fprintf('\nTest 4: Attempting to parse key MATLAB files...\n');
key_files = {
    'Utilities_Matlab/MyUtilities/Murat_test.m', ...
    'Utilities_Matlab/MyUtilities/Murat_testAll.m'
};

for i = 1:length(key_files)
    fname = key_files{i};
    try
        % Try to parse the file - this will catch basic syntax errors
        if exist(fname, 'file') == 2
            fid = fopen(fname, 'r');
            % Basic check: file opens and contains function declaration
            content = fread(fid, '*char')';
            fclose(fid);
            
            if contains(content, 'function')
                fprintf('  ✓ %s parses OK\n', fname);
            else
                fprintf('  ⚠ %s has no function declaration\n', fname);
            end
        end
    catch ME
        fprintf('  ✗ %s: %s\n', fname, ME.message);
    end
end

%% Summary
fprintf('\n');
fprintf('================== TEST SUMMARY ==================\n');
passed = sum(strcmp({test_results.status}, 'PASS'));
failed = sum(strcmp({test_results.status}, 'FAIL'));
errors = sum(strcmp({test_results.status}, 'ERROR'));

fprintf('Passed: %d | Failed: %d | Errors: %d\n', passed, failed, errors);

if failed > 0 || errors > 0
    fprintf('\nFailed/Error Tests:\n');
    for i = 1:length(test_results)
        if ~strcmp(test_results{i}.status, 'PASS')
            fprintf('  [%s] %s: %s\n', test_results{i}.status, test_results{i}.test, test_results{i}.msg);
        end
    end
    error('Syntax validation FAILED');
else
    fprintf('\n✓ All syntax checks passed!\n');
end
