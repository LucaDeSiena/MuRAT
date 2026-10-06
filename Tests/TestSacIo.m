classdef TestSacIo < matlab.unittest.TestCase
    % TESTSACIO  SAC reader availability, sample-data and header checks.
    %
    %   Replaces the script-based test_sac_io.m. Paths are anchored to the
    %   repository root, so the tests do not depend on the current folder.
    %   The header test uses the small SAC file shipped with MatSAC
    %   (Utilities_Matlab/MatSAC/N.MYJH.Z.sac), so it needs no extra data.
    %   Dataset tests are skipped (not failed) when a sample dataset folder
    %   is not present.

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

        function sacHeaderHasExpectedStructure(tc)
            % fget_sac returns [time, data, header]; check the header fields MuRAT relies on.
            tc.assumeNotEmpty(which('fget_sac'), 'fget_sac not on path.');
            sacFile = fullfile(tc.Root, 'Utilities_Matlab', 'MatSAC', 'N.MYJH.Z.sac');
            tc.assumeTrue(isfile(sacFile), 'MatSAC sample file N.MYJH.Z.sac not found.');

            [t, data, hdr] = fget_sac(sacFile);

            tc.verifyNotEmpty(data, 'No samples were read.');
            tc.verifyEqual(numel(t), numel(data), 'Time and data vectors differ in length.');
            tc.verifyFields(hdr,         {'times','station','event','user','descrip','evsta','llnl'}, 'hdr');
            tc.verifyFields(hdr.times,   {'delta','b','e','o','a','t0','t1'},                         'hdr.times');
            tc.verifyFields(hdr.event,   {'evla','evlo','evdp'},                                      'hdr.event');
            tc.verifyFields(hdr.station, {'stla','stlo','stel','stdp'},                               'hdr.station');
            tc.verifyGreaterThan(hdr.times.delta, 0, 'Sampling interval must be positive.');
        end

        function sampleDatasetHasFiles(tc, sampleFolder)
            folder = fullfile(tc.Root, sampleFolder);
            tc.assumeTrue(isfolder(folder), [sampleFolder ' is not in the repository.']);
            listing = dir(fullfile(folder, '**', '*'));
            listing = listing(~[listing.isdir]);
            listing = listing(~startsWith({listing.name}, '.'));
            tc.verifyNotEmpty(listing, ['No data files found in ' sampleFolder '.']);
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
