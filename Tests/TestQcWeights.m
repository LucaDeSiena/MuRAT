classdef TestQcWeights < matlab.unittest.TestCase
    % TESTQCWEIGHTS  Weights, selection and declustering of inverse Qc
    % measurements are consistent with what Murat_Qc stores: |R| for
    % 'Linearized' (larger is better) and sigma for 'NonLinear' (smaller
    % is better).

    properties
        Root
    end

    methods (TestClassSetup)
        function addCodeToPath(tc)
            tc.Root = fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(tc.Root, 'Utilities_Matlab'), 'IncludingSubfolders', true));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(tc.Root, 'bin')));
            set(groot, 'DefaultFigureVisible', 'off');
            tc.addTeardown(@() set(groot, 'DefaultFigureVisible', 'on'));
            tc.addTeardown(@() close('all'));
        end
    end

    methods (Test)
        function linearizedWeightIsCorrelationSquared(tc)
            w = Murat_weightQc([0.9; 0.5; 0; NaN], 'Linearized');
            tc.verifyEqual(w, [0.81; 0.25; 0; 0], 'AbsTol', 1e-12);
        end

        function nonLinearWeightIsInverseVariance(tc)
            w = Murat_weightQc([1e-3; 2e-3; 0], 'NonLinear');
            tc.verifyEqual(w, [1e6; 2.5e5; 0], 'RelTol', 1e-12);
            tc.verifyGreaterThan(w(1), w(2));   % smaller sigma, larger weight
        end

        function unknownMethodErrors(tc)
            tc.verifyError(@() Murat_weightQc(1, 'Other'), ?MException);
        end

        function linearizedRetainKeepsGoodFits(tc)
            Qm  = [0.0040; 0.0041; 0.0039; 0.0040; 0.0041];
            R   = [0.95; 0.9; 0.8; 0.2; 0.1];
            keep = Murat_retainQc(0.5, Qm, R, 'Linearized');
            tc.verifyEqual(keep, [true; true; true; false; false]);
        end

        function nonLinearRetainKeepsSmallUncertainty(tc)
            Qm  = 0.004*ones(3,1);
            sig = [1e-4; 5e-3; 1e-2];
            keep = Murat_retainQc(6e-3, Qm, sig, 'NonLinear');
            tc.verifyEqual(keep, [true; true; false]);
        end

        function declusteringKeepsHighestQuality(tc)
            % Two events in the same sub-cell recorded by the same station.
            % Depths are positive, as in Murat.rays.locationsDeg.
            loc = [10.1 50.1 6.0 0 0 0; ...
                   10.1 50.1 6.1 0 0 0];
            loc(:,4:6) = repmat([10.5 50.5 0], 2, 1);
            x = 0:1000:2000; y = 0:1000:2000; z = 0:-1000:-2000;
            keep = Murat_declustering([10 50 0], [11 51 -10], x, y, z, ...
                [0.2; 0.9], loc, 1);
            tc.verifyEqual(keep, [false; true]);
            keep = Murat_declustering([10 50 0], [11 51 -10], x, y, z, ...
                [0.9; 0.2], loc, 1);
            tc.verifyEqual(keep, [true; false]);
        end

        function nonLinearUncertaintyGrowsWithMisfit(tc)
            % sigma must grow with the noise on the coda and must not
            % depend on the amplitude scale of the envelope. Before the
            % fix sigma decreased as the misfit grew.
            [env, c0, c1, tC, sr, cf] = TestQcWeights.syntheticCoda(0);
            [q0, s0] = Murat_Qc(cf, env, c0, c1, tC, sr, 'NonLinear');
            rng(7);
            [envN, ~] = TestQcWeights.syntheticCoda(0.3);
            [~, s1] = Murat_Qc(cf, envN, c0, c1, tC, sr, 'NonLinear');
            [~, s2] = Murat_Qc(cf, 1e6*envN, c0, c1, tC, sr, 'NonLinear');
            tc.verifyEqual(q0, 0.004, 'AbsTol', 2e-5);
            tc.verifyGreaterThan(s0, 0);
            tc.verifyGreaterThan(s1, s0);
            tc.verifyEqual(s2, s1, 'RelTol', 1e-9);
        end

        function nonLinearUncertaintyFromMisfitCurve(tc)
            q = (1e-5:1e-5:1e-2)';
            E = 10 + 1e4*abs(q - 0.004);       % sharp minimum
            sSharp = Murat_uncertaintyNonLinear(q, E, 50);
            sFlat  = Murat_uncertaintyNonLinear(q, 10 + 1e3*abs(q - 0.004), 50);
            tc.verifyGreaterThan(sFlat, sSharp);
            % perfect fit: sigma falls back to the grid step, never zero
            tc.verifyEqual(Murat_uncertaintyNonLinear(q, 1e4*abs(q-0.004), 50), ...
                1e-5, 'RelTol', 1e-9);
        end

        function nonLinearUncertaintyIgnoresNaNMisfits(tc)
            % Late coda at high frequency: the predictions underflow to 0
            % for the largest trial values and give NaN misfits. They have
            % zero probability; before the fix sigma became NaN and then
            % the grid step.
            q = (1e-5:1e-5:1e-2)';
            E = 10 + 1e3*abs(q - 0.004); E(q > 0.008) = NaN;
            tc.verifyEqual(Murat_uncertaintyNonLinear(q, E, 50), ...
                Murat_uncertaintyNonLinear(q(q <= 0.008), E(q <= 0.008), 50), ...
                'RelTol', 1e-12);
        end

        function inversionWeighsSquaredResidualsByW(tc)
            % Two measurements of one cell: the weighted least-squares
            % solution is sum(w.*Qm)/sum(w), with w = 1/sigma^2.
            tc.assumeTrue(~isempty(ver('optim')) && ...
                license('test', 'Optimization_Toolbox'), ...
                'Optimization Toolbox is required for Murat_inversionQc.');
            A   = [1; 1];
            Qm  = [0.004; 0.008];
            w   = Murat_weightQc([1e-3; 2e-3], 'NonLinear');   % 4:1
            out = fullfile(tempdir, 'TestQcWeightsLS');
            tc.addTeardown(@() delete([out '*']));
            s = Murat_inversionQc(A, Qm, 'Tikhonov', 5000, 10, w, [], ...
                0, [], 0, out);
            tc.verifyEqual(s.Qc, sum(w.*Qm)/sum(w), 'RelTol', 1e-3);
        end

        function badMeasurementBarelyMovesTheModel(tc)
            % One outlier with a poor NonLinear fit (large sigma) must have
            % little effect; before the fix it had the largest weight.
            tc.assumeTrue(~isempty(ver('optim')) && ...
                license('test', 'Optimization_Toolbox'), ...
                'Optimization Toolbox is required for Murat_inversionQc.');
            rng(3);
            nRays = 80; nCells = 20;
            A = rand(nRays, nCells) .* (rand(nRays, nCells) < 0.4);
            A = A ./ sum(A, 2);
            A(sum(A,2) == 0, 1) = 1;
            A(isnan(A)) = 0;
            A = A ./ sum(A, 2);
            mTrue = 0.005*ones(nCells, 1);
            Qm  = A*mTrue;
            sig = 1e-4*ones(nRays, 1);
            Qm(1)  = 0.05;          % outlier
            sig(1) = 1e-2;          % with a poor fit
            out = fullfile(tempdir, 'TestQcWeights');
            tc.addTeardown(@() delete([out '*']));
            w = Murat_weightQc(sig, 'NonLinear');
            s = Murat_inversionQc(A, Qm, 'Tikhonov', 5000, 10, w, [], ...
                0.01, [], 0, out);
            tc.verifyLessThan(max(abs(s.Qc - mTrue))/0.005, 0.1);
        end
    end

    methods (Static)
        function [env, c0, c1, tC, sr, cf] = syntheticCoda(noise)
            % Single-scattering coda t^-1.5*exp(-2*pi*f*t/Qc), Qc^-1 = 0.004,
            % 20 s window starting 60 s after origin, optional
            % multiplicative noise.
            sr = 20; cf = 3; tC = 60; win = 20;
            t   = (tC + 1/sr : 1/sr : tC + win)';
            env = t.^-1.5 .* exp(-2*pi*cf*0.004*t);
            env = env .* (1 + noise*randn(size(env))).^2;
            c0  = 1; c1 = numel(env);
        end
    end
end
