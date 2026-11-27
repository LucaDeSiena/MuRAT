classdef TestPeakLapseTime < matlab.unittest.TestCase
    % TESTPEAKLAPSETIME  With lapseTimeMethod = 'Peak' each frequency has
    % its own coda window. Murat_Qc and Murat_body must use the window of
    % each frequency, and Murat_kernels takes one lapse time at a time.

    properties
        Root
    end

    methods (TestClassSetup)
        function addCodeToPath(tc)
            tc.Root = fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(tc.Root, 'Utilities_Matlab'), 'IncludingSubfolders', true));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(tc.Root, 'bin')));
        end
    end

    methods (Static)
        function [sp, srate] = envelopes(cf, invQc)
            % Synthetic coda envelopes t^-1.5 exp(-2 pi f t /Qc).
            srate = 100;
            t = (1/srate:1/srate:60)';
            sp = zeros(numel(t), numel(cf));
            for k = 1:numel(cf)
                sp(:,k) = t.^-1.5 .* exp(-2*pi*cf(k)*t*invQc(k));
            end
        end
    end

    methods (Test)
        function qcUsesEachFrequencyWindow(tc)
            cf = [2 8];
            invQc = [0.004 0.002];
            [sp, srate] = TestPeakLapseTime.envelopes(cf, invQc);
            % Before 24 s frequency 2 decays ten times faster, so using
            % the window of frequency 1 (10-20 s) would give a wrong value.
            [spFast, ~] = TestPeakLapseTime.envelopes(cf, 10*invQc);
            early = 1:24*srate;
            sp(early,2) = spFast(early,2) * sp(24*srate,2)/spFast(24*srate,2);
            tCoda = [10; 25];                 % one lapse time per frequency
            cS = floor(tCoda*srate - 1);
            cE = cS + floor(10*srate - 1);
            q = Murat_Qc(cf, sp, cS, cE, tCoda, srate, 'Linearized');
            tc.verifyEqual(q(:), invQc(:), 'RelTol', 1e-2);
        end

        function scalarWindowUnchanged(tc)
            cf = [2 8];
            invQc = [0.004 0.002];
            [sp, srate] = TestPeakLapseTime.envelopes(cf, invQc);
            tCoda = 15;
            cS = floor(tCoda*srate - 1);
            cE = cS + floor(10*srate - 1);
            q1 = Murat_Qc(cf, sp, cS, cE, tCoda, srate, 'Linearized');
            q2 = Murat_Qc(cf, sp, [cS; cS], [cE; cE], [tCoda; tCoda], ...
                srate, 'Linearized');
            tc.verifyEqual(q2, q1);
        end

        function bodyUsesEachFrequencyWindow(tc)
            srate = 100;
            sp = ones(3000, 2);
            sp(1001:2000, 1) = 4;             % coda of frequency 1
            sp(2001:3000, 2) = 9;             % coda of frequency 2
            [rBC, ~] = Murat_body(1, 0.1, srate, sp, 100, ...
                [1001; 2001], [2000; 3000]);
            % body window is ones, coda energies 4 and 9 per second
            tc.verifyEqual(rBC, [1/4 1/9], 'RelTol', 0.02);
        end

        function kernelsRejectVectorLapseTime(tc)
            tc.verifyError(@() Murat_kernels([10 20], [0 0 -5000], ...
                [1000 0 0], [0 0 0; 1000 1000 -1000], 3, 1, 0.5, 0.02, ...
                'Peak'), ?MException);
        end
    end
end
