function [image,mVp,xp,yp,zp,Xp,Yp,Zp]                    = ...
    Murat_image3D(X,Y,Z,V,color,sections,evestaz,x,y,z,divi,name)
% function image          = ...
%     Murat_image3D(X,Y,Z,V,color,sections,evestaz,x,y,z,name)
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
image                   =   figure('Name',name,...
    'NumberTitle','off','Position',[20,400,1200,1000],'visible','off');

stepgXYZ                =   [x(2)-x(1) y(2)-y(1) z(2)-z(1)];
divix                   =   stepgXYZ(1)/divi;
diviy                   =   stepgXYZ(2)/divi;
diviz                   =   stepgXYZ(3)/divi;

xp                      =   x(1)-divix:divix:x(end)+divix;
yp                      =   y(1)-diviy:diviy:y(end)+diviy;
zp                      =   z(1)-diviz:diviz:z(end)+diviz;
zp                      =   zp/1000;

[Xp,Yp,Zp]              =   meshgrid(xp,yp,zp);
mVp                     =   interp3(X,Y,Z,V,Xp,Yp,Zp);

slice(Xp, Yp, Zp, mVp, sections(2), sections(1), sections(3))

scale_mVp(1)            =   max(mVp(:));
scale_mVp(2)            =   abs(min(mVp(:)));
max_scale               =   max(scale_mVp);
if max_scale            <   1
    max_scale           =   round(max_scale,2);
    if max_scale        ==  0
        caxis([-1 1])
    else
        caxis([-max_scale max_scale])
    end
elseif max_scale        >   5
    max_scale           =   ceil(max_scale);
    caxis([-max_scale max_scale])
end

colormap(color);
colorbar
shading flat
hcb                     =   colorbar;
hcb.FontSize            =   18;

hold on
scatter3(evestaz(:,2),evestaz(:,1),evestaz(:,3),60,'c',...
    'MarkerEdgeColor','b', 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

scatter3(evestaz(:,5),evestaz(:,4),evestaz(:,6),60,'^',...
    'MarkerEdgeColor','m', 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

xlabel('Longitude (°)','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Latitude (°)','FontSize',16,'FontWeight','bold','Color','k')
zlabel('Altitude (km)','FontSize',16,'FontWeight','bold','Color','k')

xtick                   =   linspace(x(1),x(end),6);
ytick                   =   linspace(y(1),y(end),6);
ztick                   =   linspace(z(end),z(1),6)/1000;

xticks(xtick); set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
yticks(ytick); set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
zticks(ztick);
set(gca,'zticklabel',num2str(get(gca,'ztick')','%.2f'))
axis tight

SetFDefaults();
hold off
end
