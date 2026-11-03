classdef TestFunctional < matlab.unittest.TestCase
% TESTFUNCTIONAL  Toolbox-free functional smoke test for MuRAT4.
%
%   Verifies that the data-reading pipeline (fget_sac → Murat_test →
%   Murat_testData) runs end-to-end on the MSH sample dataset and returns
%   the correct output structure.
%
%   No MATLAB toolbox is required: Murat_test is called with
%   centralFrequencies=[] and figOutput=0, which skips every Signal
%   Processing call (butter, filtfilt, tukeywin).
%
%   The test is skipped (not failed) when the sac_MSH folder is absent,
%   so it is safe to run in environments that have only the code.

    properties
        Root        % repository root
        SacFolder   % full path to the MSH SAC dataset
    end

    methods (TestClassSetup)
        function addCodeToPath(tc)
            tc.Root = fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(tc.Root, 'Utilities_Matlab'), ...
                'IncludingSubfolders', true));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(tc.Root, 'bin')));
            tc.SacFolder = fullfile(tc.Root, 'sac_MSH');
        end
    end

    methods (Test)

        function fgetSacReadsOneSacFile(tc)
            % fget_sac must read a real SAC file and return a valid header.
            % Uses the first .SAC file found in sac_MSH.
            tc.assumeTrue(isfolder(tc.SacFolder), ...
                'sac_MSH not present — skipping.');
            files = dir(fullfile(tc.SacFolder, '*.SAC'));
            tc.assumeTrue(~isempty(files), ...
                'No .SAC files in sac_MSH — skipping.');

            sacFile = fullfile(files(1).folder, files(1).name);
            [t, data, hdr] = fget_sac(sacFile);

            tc.verifyNotEmpty(data,             'fget_sac returned empty data.');
            tc.verifyEqual(numel(t), numel(data), ...
                'Time and data vectors differ in length.');
            tc.verifyGreaterThan(hdr.times.delta, 0, ...
                'Sampling interval must be positive.');
            tc.verifyNotEqual(hdr.event.evla, -12345, ...
                'Event latitude is missing in the MSH SAC file.');
            tc.verifyNotEqual(hdr.times.a, -12345, ...
                'P-pick is missing in the MSH SAC file.');
        end

        function muratTestRunsWithoutToolboxes(tc)
            % Murat_test with centralFrequencies=[] and figOutput=0
            % must return the SAC header without calling any toolbox.
            tc.assumeTrue(isfolder(tc.SacFolder), ...
                'sac_MSH not present — skipping.');
            files = dir(fullfile(tc.SacFolder, '*.SAC'));
            tc.assumeTrue(~isempty(files), ...
                'No .SAC files in sac_MSH — skipping.');

            sacFile = fullfile(files(1).folder, files(1).name);
            % centralFrequencies=[], figOutput=0, verboseOutput=0
            [~, SAChdr] = Murat_test(sacFile, [], 8, 0, 0);

            tc.verifyClass(SAChdr, 'struct', ...
                'Murat_test did not return a struct header.');
            tc.verifyTrue(isfield(SAChdr, 'times'), ...
                'SAChdr is missing the times sub-struct.');
            tc.verifyTrue(isfield(SAChdr, 'event'), ...
                'SAChdr is missing the event sub-struct.');
            tc.verifyTrue(isfield(SAChdr, 'station'), ...
                'SAChdr is missing the station sub-struct.');
        end

        function muratTestDataBuildsTable(tc)
            % Murat_testData must process all files in sac_MSH and return
            % a table with the correct dimensions and content.
            tc.assumeTrue(isfolder(tc.SacFolder), ...
                'sac_MSH not present — skipping.');
            files = dir(fullfile(tc.SacFolder, '*.SAC'));
            tc.assumeTrue(~isempty(files), ...
                'No .SAC files in sac_MSH — skipping.');

            originTime = 'SAChdr.times.o';
            PTime      = 'SAChdr.times.a';
            STime      = 'SAChdr.times.t0';

            [header, flag, sacHeader] = ...
                Murat_testData(tc.SacFolder, originTime, PTime, STime);

            % Table structure
            tc.verifyClass(header, 'table', ...
                'Murat_testData did not return a table.');
            tc.verifyEqual(height(header), numel(files), ...
                'Table has wrong number of rows.');
            for col = {'Names','P','EvLat','EvLon','EvDepth','StLat','StLon','StElev'}
                tc.verifyTrue(any(strcmp(header.Properties.VariableNames, col{1})), ...
                    ['Table missing column: ' col{1}]);
            end

            % P-picks must all be populated (mandatory field)
            pPicks = header.P;
            tc.verifyTrue(all(~cellfun(@isempty, pPicks)), ...
                'One or more P-picks are empty — mandatory field missing.');

            % Coordinates must be populated
            tc.verifyTrue(all(~cellfun(@isempty, header.EvLat)), ...
                'One or more event latitudes are empty.');

            % flag: MSH has no origin time or S-pick → expect flag = 3
            tc.verifyEqual(flag, 3, ...
                ['flag should be 3 (origin + S missing) for MSH; got ' num2str(flag)]);

            % sacHeader struct must have one field per file
            tc.verifyClass(sacHeader, 'struct', ...
                'sacHeader output is not a struct.');
            tc.verifyEqual(numel(fieldnames(sacHeader)), numel(files), ...
                'sacHeader has wrong number of fields.');
        end

    end
end
