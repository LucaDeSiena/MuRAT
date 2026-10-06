classdef TestParallel < matlab.unittest.TestCase
    % TESTPARALLEL  parfor smoke test. Skipped if Parallel Computing Toolbox is absent.
    %   Replaces the old test.m (whose checks were trivial and asserted nothing).

    methods (Test)
        function parforComputesCorrectResult(tc)
            tc.assumeTrue(~isempty(ver('parallel')) && license('test', 'Distrib_Computing_Toolbox'), ...
                'Parallel Computing Toolbox not available.');
            n = 5;
            out = zeros(n, 1);
            parfor i = 1:n
                out(i) = i^2;
            end
            tc.verifyEqual(out, (1:n)'.^2);
        end
    end
end
