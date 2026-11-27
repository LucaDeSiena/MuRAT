classdef TestRayTracing < matlab.unittest.TestCase
    % TESTRAYTRACING  Pseudo-bending ray tracing on small synthetic grids.
    %
    %   Rays may start or end on the grid edges and may leave the grid,
    %   as for sources at the margins of a tremor-location grid.

    properties
        grid
    end

    methods (TestClassSetup)
        function addCodeToPath(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'bin')));
        end
    end

    methods (TestMethodSetup)
        function makeGrid(tc)
            tc.grid.x = 0:100:2000;
            tc.grid.y = 0:100:1500;
            tc.grid.z = 1000:-100:-2000;
        end
    end

    methods
        function pvel = gradientModel(tc)
            % pvel(y,x,z) in km/s, faster with depth
            [~, ~, Z] = meshgrid(tc.grid.x, tc.grid.y, tc.grid.z);
            pvel = 1.5 + 0.001*(1000 - Z);
        end
    end

    methods (Test)
        function velocityOnEveryGridEdge(tc)
            pvel = tc.gradientModel();
            g = tc.grid;
            corners = [g.x(1) g.y(1) g.z(1); g.x(end) g.y(end) g.z(end); ...
                g.x(1) g.y(end) g.z(1); g.x(end) g.y(1) g.z(end)];
            for k = 1:size(corners, 1)
                v = Murat_velocity(corners(k,1), corners(k,2), corners(k,3), g, pvel);
                tc.verifyEqual(v, 1.5 + 0.001*(1000 - corners(k,3)), 'AbsTol', 1e-12);
                [vx, vy, vz] = Murat_veldicity(corners(k,1), corners(k,2), corners(k,3), g, pvel);
                tc.verifyEqual([vx vy vz], [0 0 -0.001], 'AbsTol', 1e-12);
            end
        end

        function velocityContinuousOutsideGrid(tc)
            pvel = tc.gradientModel();
            g = tc.grid;
            vIn  = Murat_velocity(500, 700, -2000, g, pvel);
            vOut = Murat_velocity(500, 700, -2500, g, pvel);
            tc.verifyEqual(vOut, vIn, 'AbsTol', 1e-12);
            vSide = Murat_velocity(-300, 700, -350, g, pvel);
            tc.verifyEqual(vSide, 1.5 + 0.001*1350, 'AbsTol', 1e-12);
            [vx, vy, vz] = Murat_veldicity(500, 700, -2500, g, pvel);
            tc.verifyEqual([vx vy vz], [0 0 0], 'AbsTol', 1e-12);
        end

        function homogeneousRayIsStraight(tc)
            g = tc.grid;
            pvel = 2*ones(numel(g.y), numel(g.x), numel(g.z));
            ray = [g.x(1) g.x(end); g.y(1) g.y(end); g.z(1) g.z(end)];
            ma = Murat_tracing(ray, g, pvel);
            d = norm(ray(:,2) - ray(:,1));
            tc.verifyEqual(ma(end,5), d/2000, 'RelTol', 1e-9);
            tc.verifyEqual(ma([1 end], 2:4), ray', 'AbsTol', 1e-9);
        end

        function raysStayInsideGrid(tc)
            % Velocity grows with depth, so rays between deep points would
            %   dive below the grid bottom
            pvel = tc.gradientModel();
            g = tc.grid;
            src = [g.x(1) 700 -1900];
            rec = [g.x(end) 800 -1950; g.x(end) 800 1200];
            for k = 1:size(rec, 1)
                ma = Murat_tracing([src' rec(k,:)'], g, pvel);
                tc.verifyGreaterThanOrEqual(min(ma(:,4)), min(g.z));
                tc.verifyGreaterThanOrEqual(min(ma(:,2)), min(g.x));
                tc.verifyLessThanOrEqual(max(ma(:,2)), max(g.x));
                tc.verifyEqual(ma(end,2:4), rec(k,:), 'AbsTol', 1e-9);
            end
        end

        function edgeRaysGiveTravelTimes(tc)
            pvel = tc.gradientModel();
            g = tc.grid;
            src = [g.x(1) g.y(1) g.z(end); g.x(1) 800 g.z(end); ...
                1000 g.y(1) g.z(end); g.x(end) g.y(end) -1000];
            rec = [g.x(end) g.y(1) g.z(1); g.x(1) g.y(end) g.z(1)];
            for i = 1:size(src, 1)
                for k = 1:size(rec, 1)
                    ma = Murat_tracing([src(i,:)' rec(k,:)'], g, pvel);
                    d = norm(rec(k,:) - src(i,:));
                    tt = ma(end,5);
                    % Faster than the slowest velocity on the path, and
                    % never zero
                    tc.verifyGreaterThan(tt, d/3500);
                    tc.verifyLessThan(tt, d/1500);
                    tc.verifyEqual(ma(:,5), sort(ma(:,5)));
                end
            end
        end
    end
end
