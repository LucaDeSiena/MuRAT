function [locationsDeg_i,locationsM_i]  = Murat_location(origin,SAChdr_i,waveformName)
% function [locationsDeg_i, locationsM_i]  =   Murat_location(origin,SAChdr_i,waveformName)
%
% CREATES the variable with event and station coordinates
% CALCULATES the coordinates of event and station in meters.
% The coordinates and event info come from the SAC file.
%
% Input parameters:
%    origin:            origin of the inversion grid
%    SAChdr_i:          SAC header
%    waveformName:      name of the waveform file (optional, for error messages)
%
% Output parameters:
%    locationsDeg_i:    locations in degrees
%    locationsM_i:      locations in meters

if nargin < 3
    waveformName = 'unknown';
end

even_raw = [SAChdr_i.event.evla   SAChdr_i.event.evlo   SAChdr_i.event.evdp];
if any(even_raw == -12345)
    error(['Waveform ' waveformName ' has missing event location field!'])
end
even    =   [even_raw(1) even_raw(2) even_raw(3)*1000];

stati_raw = [SAChdr_i.station.stla SAChdr_i.station.stlo   SAChdr_i.station.stel];
if any(stati_raw == -12345)
    error(['Waveform ' waveformName ' has missing station location field!'])
end
stati   =   stati_raw;

% For plotting
locationsDeg_i      =   [even stati];

wgs84               =   wgs84Ellipsoid("m");
[distEven,azEven]   =...
                     distance(origin(1),origin(2),even(1),even(2),wgs84);
[distStati,azStati] =...
                     distance(origin(1),origin(2),stati(1),stati(2),wgs84);

% Transforms in meters
dist_x_even         =   distEven*sin(azEven*2*pi/360);
dist_y_even         =   distEven*cos(azEven*2*pi/360);
dist_x_station      =   distStati*sin(azStati*2*pi/360);
dist_y_station      =   distStati*cos(azStati*2*pi/360);
dist_z_even         =   -even(3);
dist_z_station      =   stati(3);

% For ray tracing
locationsM_i        =   [dist_x_even dist_y_even dist_z_even...
    dist_x_station dist_y_station dist_z_station];
