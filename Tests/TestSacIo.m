classdef TestSacIo < matlab.unittest.TestCase
    % TESTSACIO  SAC reader availability, sample-data and header checks.
    %
    %   Replaces the script-based test_sac_io.m. Paths are anchored to the
    %   repository root, so the tests do not depend on the current folder.
    %   Tests are skipped (not failed) when a sample dataset folder is not present.

    properties
        Root   % repository root
    end

    properties (TestParameter)
        sampleFolder = struct('MSH', 'sac_MSH', 'Romania', 'sac_Romania', 'Toba', 'sac_Toba');
    end

    methods (TestClassSetup)
        function addCodeToPath(tc)
            tc.Root = fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(tc.Root, 'Utilities_Matlab'), 'IncludingSubfolders', true));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(tc.Root, 'bin')));
        end
    end

    methods (Test)
        function sacReaderIsOnPath(tc)
            tc.verifyNotEmpty(which('fget_sac'), ...
                'fget_sac.m was not found under Utilities_Matlab.');
        end

        function sampleDatasetHasFiles(tc, sampleFolder)
            folder = fullfile(tc.Root, sampleFolder);
            tc.assumeTrue(isfolder(folder), [sampleFolder ' is not in the repository.']);
            listing = dir(fullfile(folder, '**', '*'));
            listing = listing(~[listing.isdir]);
            listing = listing(~startsWith({listing.name}, '.'));
            tc.verifyNotEmpty(listing, ['No data files found in ' sampleFolder '.']);
        end

        function sacHeaderHasExpectedStructure(tc)
            % Reads one Toba SAC file and checks the header fields MuRAT relies on.
            % Assumes the call [~,~,hdr] = fget_sac(file); check `help fget_sac`.
            tc.assumeNotEmpty(which('fget_sac'), 'fget_sac not on path.');
            folder = fullfile(tc.Root, 'sac_Toba');
            tc.assumeTrue(isfolder(folder), 'sac_Toba is not in the repository.');
            f = [dir(fullfile(folder, '**', '*.sac')); dir(fullfile(folder, '**', '*.SAC'))];
            tc.assumeNotEmpty(f, 'No .sac files found in sac_Toba.');

            [~, ~, hdr] = fget_sac(fullfile(f(1).folder, f(1).name));

            tc.verifyFields(hdr,         {'times','station','event','user','descrip','evsta','llnl'}, 'hdr');
            tc.verifyFields(hdr.times,   {'delta','b','e','o','a','t0','t1'},                         'hdr.times');
            tc.verifyFields(hdr.event,   {'evla','evlo','evdp'},                                      'hdr.event');
            tc.verifyFields(hdr.station, {'stla','stlo','stel','stdp'},                               'hdr.station');
        end
    end

    methods (Access = private)
        function verifyFields(tc, s, names, label)
            missing = names(~isfield(s, names));
            tc.verifyEmpty(missing, ...
                sprintf('%s is missing fields: %s', label, strjoin(missing, ', ')));
        end
    end
end
