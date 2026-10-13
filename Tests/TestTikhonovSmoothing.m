classdef TestTikhonovSmoothing < matlab.unittest.TestCase
    % TESTTIKHONOVSMOOTHING  First-difference smoothing in the Tikhonov
    % inversions of Q and Qc (Murat_gradientOperator, Murat_regulariser,
    % Murat_tikhonovQ, Murat_tikhonovQc, Murat_inversionQ).

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
        function [coords, dims] = grid(nx, ny, nz)
            [Z, Y, X] = ndgrid(-(0:nz-1)*2000, (0:ny-1)*8000, (0:nx-1)*9000);
            coords = [X(:) Y(:) Z(:)];
            dims = [nx ny nz];
        end

        function [A, mTrue, nodes] = smoothProblem(seed)
            % Rays as random straight segments through a 12 x 12 x 1 grid.
            rng(seed);
            nx = 12; ny = 12;
            [cx, cy] = ndgrid(1:nx, 1:ny);
            nodes = numel(cx);
            mTrue = 0.004 + 0.002 * sin(cx(:)/3) .* cos(cy(:)/4);
            nRays = 300;
            A = zeros(nRays, nodes);
            for r = 1:nRays
                p = [rand*nx rand*ny]; q = [rand*nx rand*ny];
                s = linspace(0, 1, 200)';
                xy = p + s .* (q - p);
                idx = sub2ind([nx ny], min(nx, floor(xy(:,1)) + 1), min(ny, floor(xy(:,2)) + 1));
                A(r,:) = accumarray(idx, norm(q - p) / numel(s), [nodes 1])';
            end
        end
    end

    methods (Test)
        function operatorHasOneRowPerAdjacentPair(tc)
            [coords, d] = TestTikhonovSmoothing.grid(4, 3, 5);
            perm = randperm(size(coords,1));
            D = Murat_gradientOperator(coords(perm,:));
            nPairs = (d(1)-1)*d(2)*d(3) + d(1)*(d(2)-1)*d(3) + d(1)*d(2)*(d(3)-1);
            tc.verifySize(D, [nPairs prod(d)]);
            tc.verifyEqual(full(sum(D,2)), zeros(nPairs,1));
            tc.verifyEqual(full(sum(abs(D),2)), 2*ones(nPairs,1));
            tc.verifyEqual(norm(D*ones(prod(d),1)), 0);
            c = coords(perm,:);
            [i, j] = find(D);
            rows = accumarray(i, j, [], @(v) {v});
            steps = cellfun(@(v) sum(abs(c(v(1),:) - c(v(2),:)) > 0), rows);
            tc.verifyEqual(steps, ones(nPairs,1), 'Pairs must differ along one axis only.');
        end

        function noSmoothingKeepsIdentityRegulariser(tc)
            D = Murat_gradientOperator(TestTikhonovSmoothing.grid(3, 3, 3));
            x = (1:27)';
            for sv = {[], 0}
                [L, Lt, lambda] = Murat_regulariser(27, 0.7, sv{1}, D);
                tc.verifyEqual(lambda, 0.7);
                tc.verifyEqual(L(x), x);
                tc.verifyEqual(Lt(x), x);
            end
            [L, ~, lambda] = Murat_regulariser(27, 0.7, 0.2, []);
            tc.verifyEqual(lambda, 0.7);
            tc.verifyEqual(L(x), x);
        end

        function smoothingOperatorAndAdjointAgree(tc)
            D = Murat_gradientOperator(TestTikhonovSmoothing.grid(3, 4, 2));
            n = size(D,2);
            [L, Lt, lambda] = Murat_regulariser(n, 0.3, 2, D);
            tc.verifyEqual(lambda, 1);
            x = randn(n,1); y = randn(n + size(D,1),1);
            tc.verifyEqual(L(x)'*y, x'*Lt(y), 'RelTol', 1e-12);
            tc.verifyEqual(norm(L(x))^2, 0.09*norm(x)^2 + 4*norm(D*x)^2, 'RelTol', 1e-12);
        end

        function tikhonovWithoutSmoothingIsUnchanged(tc)
            [A, mTrue] = TestTikhonovSmoothing.smoothProblem(1);
            d = A*mTrue;
            m6 = Murat_tikhonovQ(0, A, d, 0.05, mTrue, 500);
            mE = Murat_tikhonovQ(0, A, d, 0.05, mTrue, 500, [], []);
            m0 = Murat_tikhonovQc(0, A, d, 0.05, mTrue, 500, 0, ...
                Murat_gradientOperator(TestTikhonovSmoothing.grid(12, 12, 1)));
            tc.verifyEqual(mE, m6, 'AbsTol', 1e-12);
            tc.verifyEqual(m0, m6, 'AbsTol', 1e-12);
        end

        function smoothingRecoversSmoothModelBetter(tc)
            [A, mTrue] = TestTikhonovSmoothing.smoothProblem(2);
            rng(3);
            d = A*mTrue .* (1 + 0.1*randn(size(A,1),1));
            D = Murat_gradientOperator(TestTikhonovSmoothing.grid(12, 12, 1));
            lams = 10.^(-3:0.5:1);
            err = @(m) norm(m - mTrue) / norm(mTrue);
            eDamp = min(arrayfun(@(l) err(Murat_tikhonovQ(0, A, d, l, mTrue, 2000)), lams));
            eSmooth = min(arrayfun(@(l) err(Murat_tikhonovQ(0, A, d, 0, mTrue, 2000, l, D)), lams));
            tc.verifyLessThan(eSmooth, eDamp);
        end

        function smoothingKeepsMeanWhileDampingShrinksIt(tc)
            % A uniform model is in the null space of D: smoothing reduces
            % the spread but keeps the best-fitting mean, damping does not.
            [A, mTrue] = TestTikhonovSmoothing.smoothProblem(4);
            d = A*mTrue;
            D = Murat_gradientOperator(TestTikhonovSmoothing.grid(12, 12, 1));
            a1 = A*ones(size(mTrue));
            cBest = (a1'*d) / (a1'*a1);
            mWeak = Murat_tikhonovQ(0, A, d, 0, mTrue, 2000, 1, D);
            mStrong = Murat_tikhonovQ(0, A, d, 0, mTrue, 2000, 30, D);
            mDamp = Murat_tikhonovQ(0, A, d, 30, mTrue, 2000);
            tc.verifyLessThan(std(mStrong)/mean(mStrong), std(mWeak)/mean(mWeak));
            tc.verifyEqual(mean(mStrong), cBest, 'RelTol', 1e-2);
            tc.verifyLessThan(mean(mDamp), 0.5*cBest);
        end

        function inversionQUsesZeroDampingWithSmoothing(tc)
            [A, mTrue] = TestTikhonovSmoothing.smoothProblem(5);
            D = Murat_gradientOperator(TestTikhonovSmoothing.grid(12, 12, 1));
            d0 = -A*mTrue;
            Q0 = [0; mean(mTrue)];
            out = fullfile(tempdir, 'TestTikhonovSmoothing_Q');
            tc.addTeardown(@() delete([out '*']));
            [~,~,~,~,dv] = Murat_inversionQ(-A, d0, ones(size(d0)), 'Tikhonov', Q0, ...
                200, 10, [], [], 0.5, 0, out, D);
            tc.verifyEqual(dv, 0);
            [~,~,~,~,dv] = Murat_inversionQ(-A, d0, ones(size(d0)), 'Tikhonov', Q0, ...
                200, 10, [], [], [], 0, out, D);
            tc.verifyEqual(dv, 0.1);
        end
    end
end
