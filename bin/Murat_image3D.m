function [image,mVp,xp,yp,zp,Xp,Yp,Zp]                    = ...
    Murat_image3D(X,Y,Z,V,color,sections,evestaz,x,y,z,divi,name)
% function image          = ...
%     Murat_image3D(X,Y,Z,V,color,sections,evestaz,x,y,z,divi,name)
%
% PLOTS a 3D image of a field on slices.
%
% Input parameters:
%    X:         3D x matrix
%    Y:         3D y matrix
%    Z:         3D z matrix
%    V:         3D field matrix
%    color:     name of the colormap
%    sections:  location of sections
%    evestaz:   locations of earthquakes and stations in meters
%    x:         x vector
%    y:         y vector
%    z:         z vector
%    divi:      division for interpolation
%    name:      name of the title of the figure
%
% Output parameters:
%    image:     image produced
%    mVp:       interpolated field for parameter map
%    xp,yp,zp:  interpolated coords for parameter map
%    Xp,Yp,Zp:  interpolated matrices for parameter map

close all
% Create invisible figure
image                   =   myfig(name);

% Grid step and refined spacing
stepgXYZ                =   [diff(x(1:2)), diff(y(1:2)), diff(z(1:2))];
divix                   =   stepgXYZ(1)/divi;
diviy                   =   stepgXYZ(2)/divi;
diviz                   =   stepgXYZ(3)/divi;

% Build evenly spaced query vectors including one extra step on each end
nxp                     =   round((x(end)-x(1))/divix) + 1;
nyp                     =   round((y(end)-y(1))/diviy) + 1;
nzp                     =   round((z(end)-z(1))/diviz) + 1;

xp                      =   linspace(x(1)-divix, x(end)+divix, nxp);
yp                      =   linspace(y(1)-diviy, y(end)+diviy, nyp);
zp                      =   linspace(z(1)-diviz, z(end)+diviz, nzp);
zp_km                   =   zp / 1000;         % convert to km for plotting

% Create meshgrid for interpolation
[Xp,Yp,Zp]              =   meshgrid(xp,yp,zp_km);

% Interpolate (explicit method and fill missing with NaN)
% Use 'linear' by default; change to 'natural' or 'cubic' if desired
mVp = interp3(X, Y, Z, V, Xp, Yp, Zp, 'linear', NaN);

% Plot slices (note: sections is [xslice, yslice, zslice] ordering used in original)
slice(Xp, Yp, Zp, mVp, sections(2), sections(1), sections(3)); % Z in km for plotting

% Colormap and colorbar (single)
colormap(color);
hcb = colorbar;
hcb.FontSize = 18;

shading flat;
hold on;

% Plot events (earthquakes and stations)
scatter3(evestaz(:,2),evestaz(:,1),evestaz(:,3),60,'c',...
    'MarkerEdgeColor','b', 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

scatter3(evestaz(:,5),evestaz(:,4),evestaz(:,6),60,'^',...
    'MarkerEdgeColor','m', 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

% Labels and ticks
xlabel('Longitude (°)','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Latitude (°)','FontSize',16,'FontWeight','bold','Color','k')
zlabel('Altitude (km)','FontSize',16,'FontWeight','bold','Color','k')

xtick                   =   linspace(x(1),x(end),6);
ytick                   =   linspace(y(1),y(end),6);
ztick                   =   linspace(z(end),z(1),6)/1000;

xticks(xtick);  xticklabels(sprintfc('%.2f', get(gca,'xtick')'));
yticks(ytick);  yticklabels(sprintfc('%.2f', get(gca,'ytick')'));
zticks(ztick);  zticklabels(sprintfc('%.2f', get(gca,'ztick')'));

axis tight;
set(gca, 'FontSize', 14);

SetFDefaults();
hold off;
end
