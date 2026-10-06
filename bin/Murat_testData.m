function [muratHeader,flag,sacHeader] =...
    Murat_testData(folderPath,originTime,PTime,STime)
% TEST all seismograms in a folder for the input parameters and
% CREATES a file storing the parameters and flagging those missing
%
%	Input Parameters:
%       folderPath:     folder containing the SAC data
%       originTime:     origin time variable selected by the user
%       PTime:          P time variable selected by the user
%       STime:          S time variable selected by the user
%
%   Output:
%       SACheader:      SAC headers for all files;
%       muratHeader:	Murat table showing the necessary parameter
%       flag:           flags missing optional variables
%

Names               =	createsList(folderPath);
nFiles              =   numel(Names);

Origin              =   cell(nFiles,1);
P                   =   cell(nFiles,1);
S                   =   cell(nFiles,1);
EvLat               =   cell(nFiles,1);
EvLon               =   cell(nFiles,1);
EvDepth             =   cell(nFiles,1);
StLat               =   cell(nFiles,1);
StLon               =   cell(nFiles,1);
StElev              =   cell(nFiles,1);

origMissing         =   false;
sMissing            =   false;
sacHeader           =   struct();

for i=1:nFiles
    fname           =   Names{i};    
    [~,SAChdr]      =   Murat_test(fname,[],8,0,0);
    fld             =   sprintf('SAC_%g', i);
    sacHeader.(fld) =   SAChdr;

    originVal       =   getFieldValue(SAChdr, originTime);
    if isequal(eval(originTime),-12345) 
        Origin{i}   =   [];
        origMissing =   true;
    else
        Origin{i}   =   originVal;
    end
    
    pVal            =   getFieldValue(SAChdr, PTime);
    if isequal(eval(PTime),-12345)
        error(['There is a missing P value, check ' fname '!']);
    else
        P{i}        =   pVal;
    end
    
    sVal            =   getFieldValue(SAChdr, STime);
    if isequal(eval(STime),-12345)
        S{i}        =   []; 
        sMissing    =   true;
    else
        S{i}        =   sVal;
    end
    
    if isequal(SAChdr.event.evla,-12345)
        error(['There is a missing event coordinate, check ' fname '!']);
    else
        EvLat{i}    =   SAChdr.event.evla;
    end
    
    if isequal(SAChdr.event.evlo,-12345)
        error(['There is a missing event coordinate, check ' fname '!']);
    else
        EvLon{i}    =   SAChdr.event.evlo;
    end
    
    if isequal(SAChdr.event.evdp,-12345)
        error(['There is a missing event coordinate, check ' fname '!']);
    else
        EvDepth{i}  =   SAChdr.event.evdp;
    end
    
    if isequal(SAChdr.station.stla, -12345)
        error(['There is a missing station coordinate, check ' fname '!']);
    else
        StLat{i}    =   SAChdr.station.stla;
    end
    
    if isequal(SAChdr.station.stlo,-12345)
        error(['There is a missing station coordinate, check ' fname '!']);
    else
        StLon{i}    =   SAChdr.station.stlo;
    end
    
    if isequal(SAChdr.station.stel,-12345)
        error(['There is a missing station coordinate, check ' fname '!']);
    else
        StElev{i}   =   SAChdr.station.stel;
    end
    
end


muratHeader         =   table(Names,Origin,P,S,EvLat,EvLon,...
    EvDepth,StLat,StLon,StElev);

flag                =   origMissing + 2*sMissing;

end

%% Helper: Safe nested field read from SAChdr using dot-separated path
function val = getFieldValue(SAChdr, pathStr)
% pathStr examples: 'SAChdr.times.o' or 'SAChdr.times.t0'
% We ignore leading "SAChdr." if present.
if startsWith(pathStr, 'SAChdr.')
    pathStr         =   pathStr(8:end);
end
parts               =   strsplit(pathStr, '.');
val                 =   SAChdr;
for k = 1:numel(parts)
    fld             =   parts{k};
    if isstruct(val) && isfield(val, fld)
        val         =   val.(fld);
    else
        val         =   []; % missing field -> return empty (not -12345)
        return
    end
end
end

%% Helper: createsList
function [listWithFolder, listNoFolder] = createsList(directory)
d                   =   dir(directory);
if isempty(d)
    listWithFolder  =   {};
    listNoFolder    =   {};
    return
end
names               =   {d.name}.';
folders             =   {d.folder}.';
mask                =   ~startsWith(names, '.');
names               =   names(mask);
folders             =   folders(mask);
listWithFolder      =   fullfile(folders, names);
listNoFolder        =   names;
end