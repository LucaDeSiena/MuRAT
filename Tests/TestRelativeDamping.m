classdef TestRelativeDamping < matlab.unittest.TestCase
    % TESTRELATIVEDAMPING  Damping relative to norm(A) in Murat_inversionQ
    % and Murat_inversionQc, and the warm start and iteration-limit
    % warning of Murat_tikhonov_nonneg_fista.

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

    methods (Static)
        function [A, mTrue] = problem(seed)
            % Random nonnegative ray lengths through 40 cells.
            rng(seed);
            nRays = 120; nCells = 40;
            A = rand(nRays, nCells) .* (rand(nRays, nCells) < 0.3);
            mTrue = 0.004 + 0.002 * rand(nCells, 1);
        end

        function out = figureBase()
            out = fullfile(tempdir, 'TestRelativeDamping');
        end
    end

    methods (Test)
        function fistaWarnsAtIterationLimit(tc)
            state = warning('on', 'MuRAT:fistaMaxIterations');
            tc.addTeardown(@() warning(state));
            [A, mTrue] = TestRelativeDamping.problem(1);
            opts = struct('maxit', 2, 'tol', 1e-12);
            tc.verifyWarning(@() Murat_tikhonov_nonneg_fista(@(x) A*x, ...
                @(y) A'*y, @(x) x, @(y) y, A*mTrue, size(A,2), 1e-6, opts), ...
                'MuRAT:fistaMaxIterations');
        end

        function fistaStartsFromX0(tc)
            % Starting at the exact solution of a consistent problem with no
            % damping, FISTA stops at once and returns that solution.
            [A, mTrue] = TestRelativeDamping.problem(2);
            opts = struct('maxit', 5, 'tol', 1e-7, 'x0', mTrue);
            x = tc.verifyWarningFree(@() Murat_tikhonov_nonneg_fista( ...
                @(x) A*x, @(y) A'*y, @(x) x, @(y) y, A*mTrue, ...
                size(A,2), 0, opts));
            tc.verifyEqual(x, mTrue, 'RelTol', 1e-6);
        end

        function qDampingIsScaleInvariant(tc)
            tc.assumeTrue(~isempty(ver('optim')) && ...
                license('test', 'Optimization_Toolbox'), ...
                'Optimization Toolbox is required for Murat_inversionQ.');
            [A, mTrue] = TestRelativeDamping.problem(3);
            d0 = -A*mTrue;
            Q0 = [0; mean(mTrue)];
            out = TestRelativeDamping.figureBase();
            tc.addTeardown(@() delete([out '*']));
            for c = [1e-3 1e3]
                s1 = Murat_inversionQ(-A, d0, [], 'Tikhonov', Q0, 5000, 10, ...
                    [], 0.1, [], 0, out);
                s2 = Murat_inversionQ(-c*A, c*d0, [], 'Tikhonov', Q0, 5000, 10, ...
                    [], 0.1, [], 0, out);
                tc.verifyEqual(s2.Q, s1.Q, 'RelTol', 1e-4, ...
                    sprintf('Scaling A and d by %g changed the Q model.', c));
            end
        end

        function qcDampingIsScaleInvariant(tc)
            tc.assumeTrue(~isempty(ver('optim')) && ...
                license('test', 'Optimization_Toolbox'), ...
                'Optimization Toolbox is required for Murat_inversionQc.');
            [A, mTrue] = TestRelativeDamping.problem(4);
            Qm = A*mTrue ./ sum(A,2);
            RZZ = ones(size(Qm));
            out = TestRelativeDamping.figureBase();
            tc.addTeardown(@() delete([out '*']));
            s1 = Murat_inversionQc(A, Qm, 'Tikhonov', 5000, 10, RZZ, [], ...
                0.1, [], 0, out);
            % Larger data weights scale the weighted matrix and data alike.
            s2 = Murat_inversionQc(A, Qm, 'Tikhonov', 5000, 10, RZZ/30, [], ...
                0.1, [], 0, out);
            tc.verifyEqual(s2.Qc, s1.Qc, 'RelTol', 1e-4);
        end

        function moreDampingGivesSmallerModelAndLargerMisfit(tc)
            tc.assumeTrue(~isempty(ver('optim')) && ...
                license('test', 'Optimization_Toolbox'), ...
                'Optimization Toolbox is required for Murat_inversionQ.');
            [A, mTrue] = TestRelativeDamping.problem(5);
            rng(6);
            d0 = -A*mTrue .* (1 + 0.05*randn(size(A,1),1));
            Q0 = [0; mean(mTrue)];
            out = TestRelativeDamping.figureBase();
            tc.addTeardown(@() delete([out '*']));
            damps = [0.03 0.1 0.3 1];
            nrm = zeros(size(damps)); mis = zeros(size(damps));
            for i = 1:numel(damps)
                s = Murat_inversionQ(-A, d0, [], 'Tikhonov', Q0, 5000, 10, ...
                    [], damps(i), [], 0, out);
                nrm(i) = norm(s.Q);
                mis(i) = norm(-A*s.Q - d0);
            end
            tc.verifyTrue(all(diff(nrm) < 0), 'Model norm must fall as damping grows.');
            tc.verifyTrue(all(diff(mis) > 0), 'Misfit must grow with damping.');
        end
    end
end
