classdef TestUnfoldXYZ < matlab.unittest.TestCase
    % TESTUNFOLDXYZ  Murat_unfoldXYZ ordering, including single-node axes.
    %
    %   Murat_kernels can build a grid with one node along an axis when the
    %   lapse time is short compared with the model grid step.

    properties (TestParameter)
        axes3 = struct( ...
            'regular', {{[1;2;3], [10;20], [-1;-2;-3;-4]}}, ...
            'singleX', {{5, [10;20], [-1;-2]}}, ...
            'singleY', {{[1;2], 72, [-1;-2;-3]}}, ...
            'singleZ', {{[1;2], [10;20;30], -7}}, ...
            'allSingle', {{1, 2, 3}});
    end

    methods (TestClassSetup)
        function addCodeToPath(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'bin')));
        end
    end

    methods (Test)
        function matchesNdgridOrder(tc, axes3)
            [x, y, z] = axes3{:};
            r = Murat_unfoldXYZ(x, y, z);
            [Z, Y, X] = ndgrid(z, y, x);
            tc.verifyEqual(r, [X(:) Y(:) Z(:)]);
        end
    end
end
