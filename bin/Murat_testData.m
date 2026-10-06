function [muratHeader, flag, sacHeader] = ...
    Murat_testData(folderPath, originTime, PTime, STime)
% TEST all seismograms in a folder for the input parameters and
% CREATES a file storing the parameters and flagging those missing.
%
%   Input Parameters:
%       folderPath:     folder containing the SAC data
%       originTime:     dot-path to the origin-time field,
%                       e.g. 'SAChdr.times.o'
%       PTime:          dot-path to the P-pick field,
%                       e.g. 'SAChdr.times.a'
%       STime:          dot-path to the S-pick field,
%                       e.g. 'SAChdr.times.t0'
%
%   Output:
%       muratHeader:    table with one row per SAC file and columns
%                       Names, Origin, P, S, EvLat, EvLon, EvDepth,
%                       StLat, StLon, StElev
%       flag:           integer encoding of missing optional fields:
%                         0 – all present
%                         1 – origin times missing
%                         2 – S-wave times missing
%                         3 – both missing
%       sacHeader:      struct with one field per file (SAC_1, SAC_2 …)
%                       containing the full SAC header struct returned
%                       by Murat_test
%
%   MuRAT uses -12345 as the SAC missing-value sentinel.
%   Origin time and S-pick are optional (flagged); P-pick and all
%   location fields are mandatory (error on missing).

Names           =   createsList(folderPath);
nFiles          =   numel(Names);

Origin          =   cell(nFiles, 1);
P               =   cell(nFiles, 1);
S               =   cell(nFiles, 1);
EvLat           =   cell(nFiles, 1);
EvLon           =   cell(nFiles, 1);
EvDepth         =   cell(nFiles, 1);
StLat           =   cell(nFiles, 1);
StLon           =   cell(nFiles, 1);
StElev          =   cell(nFiles, 1);

origMissing     =   false;
sMissing        =   false;
sacHeader       =   struct();

for i = 1:nFiles
    fname               =   Names{i};
    [~, SAChdr]         =   Murat_test(fname, [], 8, 0, 0);
    fld                 =   sprintf('SAC_%d', i);
    sacHeader.(fld)     =   SAChdr;

    %% Origin time (optional)
    originVal           =   getFieldValue(SAChdr, originTime);
    if isequal(originVal, -12345)
        Origin{i}       =   [];
        origMissing     =   true;
    else
        Origin{i}       =   originVal;
    end

    %% P-wave pick (mandatory)
    pVal                =   getFieldValue(SAChdr, PTime);
    if isequal(pVal, -12345)
        error('MuRAT:missingPick', ...
            'Missing P value in file:\n  %s', fname);
    end
    P{i}                =   pVal;

    %% S-wave pick (optional)
    sVal                =   getFieldValue(SAChdr, STime);
    if isequal(sVal, -12345)
        S{i}            =   [];
        sMissing        =   true;
    else
        S{i}            =   sVal;
    end

    %% Event location (mandatory)
    if isequal(SAChdr.event.evla, -12345)
        error('MuRAT:missingCoordinate', ...
            'Missing event latitude (evla) in file:\n  %s', fname);
    end
    EvLat{i}            =   SAChdr.event.evla;

    if isequal(SAChdr.event.evlo, -12345)
        error('MuRAT:missingCoordinate', ...
            'Missing event longitude (evlo) in file:\n  %s', fname);
    end
    EvLon{i}            =   SAChdr.event.evlo;

    if isequal(SAChdr.event.evdp, -12345)
        error('MuRAT:missingCoordinate', ...
            'Missing event depth (evdp) in file:\n  %s', fname);
    end
    EvDepth{i}          =   SAChdr.event.evdp;

    %% Station location (mandatory)
    if isequal(SAChdr.station.stla, -12345)
        error('MuRAT:missingCoordinate', ...
            'Missing station latitude (stla) in file:\n  %s', fname);
    end
    StLat{i}            =   SAChdr.station.stla;

    if isequal(SAChdr.station.stlo, -12345)
        error('MuRAT:missingCoordinate', ...
            'Missing station longitude (stlo) in file:\n  %s', fname);
    end
    StLon{i}            =   SAChdr.station.stlo;

    if isequal(SAChdr.station.stel, -12345)
        error('MuRAT:missingCoordinate', ...
            'Missing station elevation (stel) in file:\n  %s', fname);
    end
    StElev{i}           =   SAChdr.station.stel;

end

muratHeader     =   table(Names, Origin, P, S, EvLat, EvLon, ...
                          EvDepth, StLat, StLon, StElev);

% flag encodes which optional fields are absent:
%   bit 0 (value 1): origin times missing
%   bit 1 (value 2): S-wave times missing
flag            =   origMissing + 2 * sMissing;

end

%% -------------------------------------------------------------------------
%  Helper: safe nested field read from SAChdr using a dot-separated path
%  -------------------------------------------------------------------------
function val = getFieldValue(SAChdr, pathStr)
% GETFIELDVALUE  Read a field from SAChdr by dot-path string.
%
%   val = getFieldValue(SAChdr, 'SAChdr.times.o')  or
%   val = getFieldValue(SAChdr, 'times.o')
%
%   Returns the field value, or [] if any segment is absent.
%   Never calls eval().

if startsWith(pathStr, 'SAChdr.')
    pathStr     =   pathStr(8:end);      % strip leading 'SAChdr.'
end
parts           =   strsplit(pathStr, '.');
val             =   SAChdr;
for k = 1:numel(parts)
    fld         =   parts{k};
    if isstruct(val) && isfield(val, fld)
        val     =   val.(fld);
    else
        val     =   [];                  % field absent -> return empty
        return
    end
end
end

%% -------------------------------------------------------------------------
%  Helper: list visible files in a directory
%  -------------------------------------------------------------------------
function [listWithFolder, listNoFolder] = createsList(directory)
% CREATESLIST  Return full and bare names of visible files in directory.

d               =   dir(directory);
if isempty(d)
    listWithFolder  =   {};
    listNoFolder    =   {};
    return
end
names           =   {d.name}.';
folders         =   {d.folder}.';
mask            =   ~startsWith(names, '.');
names           =   names(mask);
folders         =   folders(mask);
listWithFolder  =   fullfile(folders, names);
listNoFolder    =   names;
end
