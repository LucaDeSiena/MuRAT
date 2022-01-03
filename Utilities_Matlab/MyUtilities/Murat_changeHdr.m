%% CHANGES the header of sac files to include the pickings for MSH
function seism                  =   Murat_changeHdr(newFolder)
% function seism                  =   Murat_changeHdr(newfolder)
% CHANGES header of a file
%
%	Input Parameters:
%       newFolder:              folder where you save the changed file
%
%   Output:
%       seism:                  new seismogram
%


[file,path]                     =   uigetfile('*.*');

if isequal(file,0)
   error('User selected Cancel!');   
else
   disp(['User selected ', fullfile(path,file)]);
end

seism                           =   fread_sac(fullfile(path,file));
hs                              =   seism.hdr;

Field                           =...
    ["Event Latitude";"Event Longitude";"Event depth (km)";...
    "Station Latitude";"Station Longitude"; "Station elevation (m)";...
    "Origin time"; "P time"; "S time"];

Value                           =...
    [hs.evla;hs.evlo;hs.evdp;...
    hs.stla;hs.stlo;hs.stel;...
    hs.o;hs.a;hs.t(1)];

requiredFields                  =   table(Field,Value);
disp(requiredFields);

changeAsk                       =...
    input('What do you want to change?\n Nothing?\n Event (1)?\n Station (2)?\n Time (3)?');

switch changeAsk
    case 1
        hs.evla                 =   input('Event latitude?');
        hs.evlo                 =   input('Event longitude?');
        hs.evdp                 =   input('Event depth?');
    
    case 2
        hs.stla                 =   input('Station latitude?');
        hs.stlo                 =   input('Station longitude?');
        hs.stel                 =   input('Station elevation?');

    case 3
        hs.o                    =   input('Origin time?');
        hs.a                    =   input('P time?');
        hs.t(1)                 =   input('S time?');
end

seism.hdr                       =   hs;
fwrite_sac(seism,fullfile(newFolder,file));