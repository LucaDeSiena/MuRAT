classdef TestCodaCheck < matlab.unittest.TestCase
    % TESTCODACHECK  Coda window position for every lapseTimeMethod.
    %
    %   Synthetic trace starting at t = 0 s with a non-zero origin time, so
    %   that adding the origin time twice would move the window.

    properties (TestParameter)
        originTime = {0, 30};
    end

    properties (Constant)
        Srate   = 50;
        Pick    = 8;     % S travel time from the origin (s)
        TCm     = 1.5;
        TWm     = 30;
        PD      = 2;     % peak delay (s)
    end

    methods (TestClassSetup)
        function addCodeToPath(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'bin')));
        end
    end

    methods (Test)
        function travelStartsAtMultipleOfTravelTime(tc, originTime)
            tc.verifyWindow(originTime, 'Travel', tc.TCm*tc.Pick);
        end

        function constantStartsAtFixedLapseTime(tc, originTime)
            tc.verifyWindow(originTime, 'Constant', tc.TCm);
        end

        function peakStartsAfterPeakDelay(tc, originTime)
            tc.verifyWindow(originTime, 'Peak', tc.Pick + tc.PD);
        end
    end

    methods (Access = private)
        function verifyWindow(tc, originTime, method, expectedLapse)
            tempis = (0:1/tc.Srate:180)';
            [tCoda, c0, c1] = Murat_codaCheck(originTime, originTime + tc.Pick, ...
                tc.Srate, tc.TCm, tc.TWm, tempis, tc.PD, method);
            tc.verifyEqual(tCoda, expectedLapse, 'AbsTol', 1e-12, ...
                'tCoda should be the lapse time from the origin.');
            tc.verifyEqual(tempis(c0), originTime + expectedLapse, 'AbsTol', 3/tc.Srate, ...
                'Coda window starts at the wrong sample.');
            tc.verifyEqual(tempis(c1) - tempis(c0), tc.TWm, 'AbsTol', 2/tc.Srate);
        end
    end
end
