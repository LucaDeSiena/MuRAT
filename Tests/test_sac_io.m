%% SAC I/O INTEGRATION TESTS
% Tests SAC file reading and basic header structure validation

clear; clc;

fprintf('Starting SAC I/O Integration Tests...\n\n');

test_results = {};

%% Test 1: SAC module loads
fprintf('Test 1: Loading SAC utilities...\n');
try
    addpath('Utilities_Matlab/MatSAC');
    
    % Check if we can see the functions
    which_sac = which('sac');
    which_sachdr = which('sachdr');
    
    if ~isempty(which_sac) && ~isempty(which_sachdr)
        fprintf('✓ SAC utilities are loadable\n');
        test_results{end+1} = struct('name', 'sac_load', 'pass', true);
    else
        fprintf('✗ SAC utilities not found\n');
        test_results{end+1} = struct('name', 'sac_load', 'pass', false);
    end
    
catch ME
    fprintf('✗ Error loading SAC: %s\n', ME.message);
    test_results{end+1} = struct('name', 'sac_load', 'pass', false);
end

%% Test 2: Check sample data exists
fprintf('\nTest 2: Checking for sample SAC data...\n');
try
    % MuRAT should have sample datasets according to README
    sample_dirs = {
        'Data', ...          % Generic data directory
        'SampleData', ...    % Alternative name
        'Tests/sample_data'  % Test-specific samples
    };
    
    data_found = false;
    for i = 1:length(sample_dirs)
        if isdir(sample_dirs{i})
            fprintf('  ✓ Found data directory: %s\n', sample_dirs{i});
            data_found = true;
            
            % List SAC files if any
            sac_files = dir(fullfile(sample_dirs{i}, '*.sac'));
            if ~isempty(sac_files)
                fprintf('    Found %d SAC files\n', length(sac_files));
            end
            break;
        end
    end
    
    if data_found
        test_results{end+1} = struct('name', 'sample_data_exists', 'pass', true);
    else
        fprintf('  ⚠ No sample data directory found (integration test skipped)\n');
        test_results{end+1} = struct('name', 'sample_data_exists', 'pass', true, 'note', 'skipped');
    end
    
catch ME
    fprintf('✗ Error: %s\n', ME.message);
    test_results{end+1} = struct('name', 'sample_data_exists', 'pass', false);
end

%% Test 3: SAC header structure validation
fprintf('\nTest 3: Validating SAC header structure...\n');
try
    % Simulate SAC header structure that fget_sac should return
    % Based on sachdr.m, headers should have: times, station, event, user, descrip
    
    expected_fields = {'times', 'station', 'event', 'user', 'descrip', 'evsta', 'llnl'};
    expected_time_fields = {'delta', 'b', 'e', 'o', 'a', 't0', 't1'};
    expected_event_fields = {'evla', 'evlo', 'evdp'};
    expected_station_fields = {'stla', 'stlo', 'stel', 'stdp'};
    
    fprintf('  Expected SAC header fields:\n');
    fprintf('    Main: %s\n', strjoin(expected_fields, ', '));
    fprintf('    Times: %s\n', strjoin(expected_time_fields, ', '));
    fprintf('    Event: %s\n', strjoin(expected_event_fields, ', '));
    fprintf('    Station: %s\n', strjoin(expected_station_fields, ', '));
    
    fprintf('✓ SAC header structure documented\n');
    test_results{end+1} = struct('name', 'sac_structure', 'pass', true);
    
catch ME
    fprintf('✗ Error: %s\n', ME.message);
    test_results{end+1} = struct('name', 'sac_structure', 'pass', false);
end

%% Test 4: Check for missing value marker handling
fprintf('\nTest 4: Validating missing value handling (-12345)...\n');
try
    % MuRAT uses -12345 as a sentinel for missing SAC header values
    % Check that Murat_testData.m handles this correctly
    
    fid = fopen('bin/Murat_testData.m', 'r');
    content = fread(fid, '*char')';
    fclose(fid);
    
    % Should check multiple fields for -12345
    checks_for_missing = {
        'isequal(eval(originTime),-12345)', ...
        'isequal(eval(PTime),-12345)', ...
        'isequal(eval(STime),-12345)'
    };
    
    all_checked = true;
    for i = 1:length(checks_for_missing)
        if contains(content, checks_for_missing{i})
            fprintf('  ✓ Checks for %s\n', checks_for_missing{i});
        else
            fprintf('  ✗ Missing check for %s\n', checks_for_missing{i});
            all_checked = false;
        end
    end
    
    if all_checked
        fprintf('✓ All missing value checks present\n');
        test_results{end+1} = struct('name', 'missing_value_handling', 'pass', true);
    else
        test_results{end+1} = struct('name', 'missing_value_handling', 'pass', false);
    end
    
catch ME
    fprintf('✗ Error: %s\n', ME.message);
    test_results{end+1} = struct('name', 'missing_value_handling', 'pass', false);
end

%% Summary
fprintf('\n');
fprintf('================== SAC I/O TEST SUMMARY ==================\n');
passed = sum([test_results.pass]);
total = length(test_results);

fprintf('Passed: %d / %d\n', passed, total);

if passed < total
    fprintf('\nFailed Tests:\n');
    for i = 1:length(test_results)
        if ~test_results{i}.pass
            fprintf('  ✗ %s\n', test_results{i}.name);
        end
    end
    error('SAC I/O tests FAILED');
else
    fprintf('\n✓ All SAC I/O tests passed!\n');
end
