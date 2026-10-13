classdef TestInputMethods < matlab.unittest.TestCase
    % TESTINPUTMETHODS  Example inputs only use inversion methods the code accepts.

    properties
        Root     % repository root
        Methods  % case labels in bin/Murat_inversionQ.m
    end

    properties (TestParameter)
        inputFile = {'Murat_inputToba.m', 'Murat_inputMSH.m', 'Murat_inputRomania.m'};
    end

    methods (TestClassSetup)
        function readMethods(tc)
            tc.Root = fileparts(fileparts(mfilename('fullpath')));
            src = fileread(fullfile(tc.Root, 'bin', 'Murat_inversionQ.m'));
            tc.Methods = regexp(src, '^\s*case\s+''(\w+)''', 'tokens', 'lineanchors');
            tc.Methods = [tc.Methods{:}];
        end
    end

    methods (Test)
        function inversionQcHasSameMethods(tc)
            src = fileread(fullfile(tc.Root, 'bin', 'Murat_inversionQc.m'));
            m = regexp(src, '^\s*case\s+''(\w+)''', 'tokens', 'lineanchors');
            tc.verifyEqual(sort([m{:}]), sort(tc.Methods));
        end

        function inputUsesKnownMethod(tc, inputFile)
            src = fileread(fullfile(tc.Root, inputFile));
            m = regexp(src, 'inversionMethod\s*=\s*''(\w+)''', 'tokens', 'once');
            tc.assertNotEmpty(m, ['No inversionMethod in ' inputFile '.']);
            tc.verifyTrue(ismember(m{1}, tc.Methods), ...
                sprintf('%s sets unknown inversionMethod ''%s''.', inputFile, m{1}));
        end

        function inputDocsNameKnownMethods(tc, inputFile)
            src = fileread(fullfile(tc.Root, inputFile));
            src = extractBetween(src, '## INVERSION', 'inversionMethod');
            tc.assertNotEmpty(src, ['No INVERSION section in ' inputFile '.']);
            named = regexp(src{1}, '\(\*\*''?([A-Za-z]\w*)''?\*\*\)', 'tokens');
            named = [named{:}];
            tc.assertNotEmpty(named, ['No method names documented in ' inputFile '.']);
            unknown = setdiff(named, tc.Methods);
            tc.verifyEmpty(unknown, sprintf('%s documents unknown methods: %s', ...
                inputFile, strjoin(unknown, ', ')));
        end
    end
end
