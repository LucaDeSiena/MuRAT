classdef TestHeadlessRun < matlab.unittest.TestCase
    % TESTHEADLESSRUN  MuRAT.m with a predefined inputFile (no uigetfile),
    % and the regularisation summary printed by the inversions.

    properties
        Root
    end

    methods (TestClassSetup)
        function addCodeToPath(tc)
            tc.Root = fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(tc.Root, 'bin')));
        end
    end

    methods (Test)
        function headlessReadsExtensionFromInputFile(tc)
            % An unsupported extension must reach the extension check,
            % which needs the extension of inputFile, not of 'file'.
            inputFile = [tempname '.txt']; %#ok<NASGU> used by MuRAT.m
            fid = fopen(inputFile, 'w'); fclose(fid);
            tc.addTeardown(@() delete(inputFile));
            tc.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(tc.Root));
            try
                evalc('MuRAT');
                tc.verifyFail('MuRAT.m should reject a .txt input file.');
            catch err
                tc.verifySubstring(err.message, 'Unsupported input file extension');
            end
        end

        function printsAbsoluteWeights(tc)
            out = evalc('Murat_printRegularisation(''Q'',27.5,0,true,0.1)');
            tc.verifySubstring(out, 'norm(A) = 27.5');
            tc.verifySubstring(out, 'smoothing 0.1 (absolute 2.75)');
            out = evalc('Murat_printRegularisation(''Qc'',2,0.1,false,[])');
            tc.verifySubstring(out, 'damping 0.1 (absolute 0.2)');
            tc.verifyEmpty(strfind(out, 'smoothing'));
        end
    end
end
