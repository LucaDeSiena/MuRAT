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
end
