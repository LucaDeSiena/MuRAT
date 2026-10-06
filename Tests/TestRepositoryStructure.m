classdef TestRepositoryStructure < matlab.unittest.TestCase
    % TESTREPOSITORYSTRUCTURE  Repository-level sanity checks.
    %
    %   Replaces the script-based test_syntax_check.m. Checks that critical
    %   files exist, that all sources parse, that no editor backups or
    %   duplicate Murat_* functions are committed, and guards two known
    %   regressions in Murat_testData / Murat_testAll.

    properties
        Root   % repository root
    end

    properties (TestParameter)
        criticalFile = struct( ...
            'MuRAT',           'MuRAT.m', ...
            'Murat_test',      'bin/Murat_test.m', ...
            'Murat_testData',  'bin/Murat_testData.m', ...
            'Murat_checks',    'bin/Murat_checks.m', ...
            'Murat_inversion', 'bin/Murat_inversion.m', ...
            'Murat_plot',      'bin/Murat_plot.m', ...
            'Murat_testAll',   'Utilities_Matlab/MyUtilities/Murat_testAll.m', ...
            'fget_sac',        'Utilities_Matlab/MatSAC/fget_sac.m', ...
            'README',          'README.md', ...
            'LICENSE',         'LICENSE.md');
    end

    methods (TestClassSetup)
        function locateRoot(tc)
            tc.Root = fileparts(fileparts(mfilename('fullpath')));
        end
    end

    methods (Test)
        function criticalFileExists(tc, criticalFile)
            tc.verifyTrue(isfile(fullfile(tc.Root, criticalFile)), ...
                [criticalFile ' is missing.']);
        end

        function allSourceFilesParse(tc)
            % Static parse of repo-owned .m files (no execution, no toolboxes).
            % Uses mtree, which ships with base MATLAB; ERR nodes mark parse errors.
            files = [ ...
                dir(fullfile(tc.Root, '*.m')); ...
                dir(fullfile(tc.Root, 'bin', '*.m')); ...
                dir(fullfile(tc.Root, 'Utilities_Matlab', 'MyUtilities', '*.m'))];
            tc.assertNotEmpty(files, 'No .m files found to parse.');
            bad = {};
            for k = 1:numel(files)
                f = fullfile(files(k).folder, files(k).name);
                try
                    T = mtree(f, '-file');
                    if ~isnull(mtfind(T, 'Kind', 'ERR'))
                        bad{end+1} = files(k).name; %#ok<AGROW>
                    end
                catch ME
                    bad{end+1} = [files(k).name ' (' ME.message ')']; %#ok<AGROW>
                end
            end
            tc.verifyEmpty(bad, ['Files with parse errors: ' strjoin(bad, ', ')]);
        end

        function repositoryHasNoEditorBackupFiles(tc)
            backups = dir(fullfile(tc.Root, '**', '*~'));
            names = {backups.name};
            tc.verifyEmpty(names, ['Backup files committed: ' strjoin(names, ', ')]);
        end

        function noDuplicateMuratFunctionNames(tc)
            % A duplicated name means one file silently shadows the other on the path.
            listing = [ ...
                dir(fullfile(tc.Root, 'bin', '**', 'Murat_*.m')); ...
                dir(fullfile(tc.Root, 'Utilities_Matlab', '**', 'Murat_*.m'))];
            tc.assertNotEmpty(listing, 'No Murat_*.m functions found.');
            [u, ~, idx] = unique({listing.name});
            counts = accumarray(idx(:), 1);
            dupes = u(counts > 1);
            tc.verifyEmpty(dupes, ['Function names defined more than once: ' strjoin(dupes, ', ')]);
        end

        function testDataDeclaresAndInitialisesFlag(tc)
            src = tc.readStripped('bin/Murat_testData.m');
            tc.verifyTrue(contains(src, 'function[muratHeader,flag]'), ...
                'Murat_testData must declare the output [muratHeader, flag].');
            tc.verifyTrue(contains(src, 'flag=[]'), ...
                'Murat_testData must initialise flag.');
        end

        function testAllInitialisesFlag(tc)
            src = tc.readStripped('Utilities_Matlab/MyUtilities/Murat_testAll.m');
            tc.verifyTrue(contains(src, 'flag=[]') || contains(src, 'flag=0'), ...
                'Murat_testAll must initialise flag before the loop.');
        end

        function missingValueMarkersAreChecked(tc)
            % MuRAT uses -12345 as the missing-value marker in SAC headers.
            src = tc.readStripped('bin/Murat_testData.m');
            required = {'isequal(eval(originTime),-12345)', ...
                        'isequal(eval(PTime),-12345)', ...
                        'isequal(eval(STime),-12345)'};
            for k = 1:numel(required)
                tc.verifyTrue(contains(src, required{k}), ...
                    ['Murat_testData no longer contains the check ' required{k}]);
            end
        end
    end

    methods (Access = private)
        function src = readStripped(tc, relPath)
            p = fullfile(tc.Root, relPath);
            tc.assertTrue(isfile(p), [relPath ' is missing.']);
            src = regexprep(fileread(p), '\s', '');   % ignore whitespace differences
        end
    end
end
