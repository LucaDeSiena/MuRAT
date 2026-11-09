classdef TestEnvelope < matlab.unittest.TestCase
    % TESTENVELOPE  Envelope smoothing time is honoured end to end.
    %
    %   Uses the first SAC file of the Toba sample dataset; the numerical
    %   test is skipped when that folder is not present.

    properties
        Root   % repository root
    end

    properties (TestParameter)
        dataFile = {'Murat_data.m', 'Murat_dataParallelized.m'};
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
        function smoothTimeChangesEnvelope(tc)
            % Murat_envelope uses tukeywin; CI installs base MATLAB only.
            tc.assumeTrue(~isempty(ver('signal')) && ...
                license('test', 'Signal_Toolbox'), ...
                'Signal Processing Toolbox is required for Murat_envelope.');
            listing = dir(fullfile(tc.Root, 'sac_Toba', '*.sac'));
            tc.assumeNotEmpty(listing, 'sac_Toba sample data not found.');
            sacFile = fullfile(listing(1).folder, listing(1).name);
            cf = [3 6];

            [t, spDefault] = Murat_envelope(cf, sacFile);
            [~, spOne]     = Murat_envelope(cf, sacFile, 1);
            [~, spShort]   = Murat_envelope(cf, sacFile, 0.25);

            tc.verifySize(spDefault, [numel(t) numel(cf)]);
            tc.verifyEqual(spOne, spDefault, ...
                'An explicit 1 s window should match the default.');
            tc.verifyNotEqual(spShort, spDefault, ...
                'A 0.25 s window should change the envelope.');
        end

        function dataPassesSmoothTime(tc, dataFile)
            src = fileread(fullfile(tc.Root, 'bin', dataFile));
            tc.verifyTrue(contains(src, 'Murat_envelope(cf,listSac_i,envSmoothTime)'), ...
                [dataFile ' does not pass envelopeSmoothTime to Murat_envelope.']);
        end
    end
end
