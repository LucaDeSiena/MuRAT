function muratHeader        =   Murat_testAll(folderPath)
% TEST all seismograms in a folder for the input parameters and
% CREATES a file storing the parameters and flagging those missing
%
%	Input Parameters:
%       folderPath:     folderPath
%
%   Output:
%       muratHeader:	Murat table showing the necessary parameter
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

for i=1:lengthData
    listSac_i               =   Names{i};    
    [~,SAChdr]              =   Murat_test(listSac_i,[],8,0,0);

    if isequal(SAChdr.times.o,-12345) 
        Origin{i}           =   [];
    else
        Origin{i}           =   SAChdr.times.o;
    end
    
    if isequal(SAChdr.times.a,-12345)
        P{i}                =   []; 
    else
        P{i}                =   SAChdr.times.a;
    end
    
    if isequal(SAChdr.times.t0,-12345) 
        S{i}                =   []; 
    else
        S{i}                =   SAChdr.times.t0;
    end
    
    if isequal(SAChdr.event.evla,-12345)
        EvLat{i}            =   []; 
    else
        EvLat{i}            =   SAChdr.event.evla;
    end
    
    if isequal(SAChdr.event.evlo,-12345)
        EvLon{i}            =   []; 
    else
        EvLon{i}            =   SAChdr.event.evlo;
    end
    
    if isequal(SAChdr.event.evdp,-12345)
        EvDepth{i}          =   []; 
    else
        EvDepth{i}          =   SAChdr.event.evdp;
    end
    
    if isequal(SAChdr.station.stla,-12345)
        StLat{i}            =   []; 
    else
        StLat{i}            =   SAChdr.station.stla;
    end
    
    if isequal(SAChdr.station.stlo,-12345)
        StLon{i}            =   []; 
    else
        StLon{i}            =   SAChdr.station.stlo;
    end
    
    if isequal(SAChdr.station.stel,-12345)
        StElev{i}           =   []; 
    else
        StElev{i}           =   SAChdr.station.stel;
    end
    
end

muratHeader                 =   table(Names,Origin,P,S,EvLat,EvLon,...
    EvDepth,StLat,StLon,StElev);

writetable(muratHeader,'DataHeaders.xls') 

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

end