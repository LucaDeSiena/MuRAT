function [locationsDeg_i, locationsM_i]...
                            =   Murat_location(origin,SAChdr_i)
% function [locationsDeg_i, locationsM_i]  =   Murat_location(origin,SAChdr_i)
%
% CREATES the variable with event and station coordinates
% CALCULATES the coordinates of event and station in meters.
% The coordinates and event info come from the SAC file.
%
% Input parameters:
%    origin:            origin of the inversion grid
%    SAChdr_i:          SAC header
%    
% Output parameters:
%    locationsDeg_i:    locations in degrees
%    locationsM_i:      locations in meters


even                        =   [SAChdr_i.event.evla...
    SAChdr_i.event.evlo SAChdr_i.event.evdp*1000];

stati                       =   [SAChdr_i.station.stla...
    SAChdr_i.station.stlo SAChdr_i.station.stel];

if find(even == -12345)
    error(['Waveform ' listaSac_i 'has missing event location field!'])
end

if find(stati == -12345)
    error(['Waveform ' listaSac_i 'has missing station location field!'])
end

% For plotting
locationsDeg_i              =   [even stati];

dist_xdeg_even              =   even(2)-origin(2);
dist_ydeg_even              =   even(1)-origin(1);
dist_z_even                 =   -even(3);

dist_xdeg_station           =   stati(2)-origin(2);
dist_ydeg_station           =   stati(1)-origin(1);
dist_z_station              =   stati(3);

% Transforms in meters
dist_x_even                 =   deg2km(dist_xdeg_even)*1000;
dist_y_even                 =   deg2km(dist_ydeg_even)*1000;
dist_x_station              =   deg2km(dist_xdeg_station)*1000;
dist_y_station              =   deg2km(dist_ydeg_station)*1000;

% For ray tracing
locationsM_i                =   [dist_x_even dist_y_even dist_z_even...
    dist_x_station dist_y_station dist_z_station];
