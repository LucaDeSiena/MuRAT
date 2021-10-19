function [muratHeader,flag] =...
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
%       muratHeader:	Murat table showing the necessary parameter
%       flag:           flags missing optional variables
%

[Names,~]                   =	createsList(folderPath);
lengthData                  =   length(Names);
Origin                      =   cell(lengthData,1);
P                           =   cell(lengthData,1);
S                           =   cell(lengthData,1);
EvLat                       =   cell(lengthData,1);
EvLon                       =   cell(lengthData,1);
EvDepth                     =   cell(lengthData,1);
StLat                       =   cell(lengthData,1);
StLon                       =   cell(lengthData,1);
StElev                      =   cell(lengthData,1);
flag                        =   [];
for i=1:lengthData
    listSac_i               =   Names{i};    
    [~,SAChdr]              =   Murat_test(listSac_i,[],8,0,0);

    if isequal(eval(originTime),-12345) 
        Origin{i}           =   [];
        flag                =   1;
    else
        Origin{i}           =   eval(originTime);
    end
    
    if isequal(eval(PTime),-12345)
        P{i}                =   []; 
        warning(['There is a missing P-value, check ' listSac_i '!']);
    else
        P{i}                =   eval(PTime);
    end
    
    if isequal(eval(STime),-12345) 
        S{i}                =   []; 
        flag                =   2;
    else
        S{i}                =   eval(STime);
    end
    
    if isequal(SAChdr.event.evla,-12345)
        EvLat{i}            =   []; 
        warning(['There is a missing event coordinate, check '...
            listSac_i '!']);
    else
        EvLat{i}            =   SAChdr.event.evla;
    end
    
    if isequal(SAChdr.event.evlo,-12345)
        EvLon{i}            =   []; 
        warning(['There is a missing event coordinate, check '...
            listSac_i '!']);
    else
        EvLon{i}            =   SAChdr.event.evlo;
    end
    
    if isequal(SAChdr.event.evdp,-12345)
        EvDepth{i}          =   []; 
        warning(['There is a missing event coordinate, check '...
            listSac_i '!']);
    else
        EvDepth{i}          =   SAChdr.event.evdp;
    end
    
    if isequal(SAChdr.station.stla,-12345)
        StLat{i}            =   []; 
        warning(['There is a missing station coordinate, check '...
            listSac_i '!']);
    else
        StLat{i}            =   SAChdr.station.stla;
    end
    
    if isequal(SAChdr.station.stlo,-12345)
        StLon{i}            =   []; 
        warning(['There is a missing station coordinate, check '...
            listSac_i '!']);
    else
        StLon{i}            =   SAChdr.station.stlo;
    end
    
    if isequal(SAChdr.station.stel,-12345)
        StElev{i}           =   []; 
        warning(['There is a missing station coordinate, check '...
            listSac_i '!']);
    else
        StElev{i}           =   SAChdr.station.stel;
    end
    
end


muratHeader                 =   table(Names,Origin,P,S,EvLat,EvLon,...
    EvDepth,StLat,StLon,StElev);

end
%%
function [listWithFolder,listNoFolder]...
                            =   createsList(directory)
% CREATES a list of visible files in a folder, outputs both with and
% without folder

list                        =   dir(directory);
list                        =   list(~startsWith({list.name}, '.'));

listWithFolder              =	fullfile({list.folder},{list.name})';
listNoFolder                =   {list.name}';
end