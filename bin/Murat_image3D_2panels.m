function image          = ...
    Murat_image3D_2panels(X,Y,Z,V,color,sections,evestaz,x,y,z)
% function image          = ...
%     Murat_image3D(X,Y,Z,V,color,sections,evestaz,x,y,z)
%
% PLOTS a 3D image of a field on slices on two subpanels.
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
%
% Output parameters:
%    image:     image produced in one panel

divi        =   5;
dx          =   (x(2)-x(1))/divi;
dy          =   (y(2)-y(1))/divi;
dz          =   (z(2)-z(1))/divi;

xp          =   x(1)-dx:dx:x(end)+dx;
yp          =   y(1)-dy:dy:y(end)+dy;
zp          =   (z(1)-dz:dz:z(end)+dz)/1000;

[Xp,Yp,Zp]  =   meshgrid(xp,yp,zp);
mVp         =   interp3(X,Y,Z,V,Xp,Yp,Zp);
z           =   sort(z)/1000;

image       =slice(Xp, Yp, Zp, mVp, sections(2), sections(1), sections(3));

colormap(color);
colorbar
shading flat
hcb         =   colorbar;
hcb.FontSize=   18;

hold on
scatter3(evestaz(:,2),evestaz(:,1),evestaz(:,3),60,'c',...
    'MarkerEdgeColor','b', 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

scatter3(evestaz(:,5),evestaz(:,4),evestaz(:,6),60,'^',...
    'MarkerEdgeColor','m', 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

xlabel('Longitude (°)','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Latitude (°)','FontSize',16,'FontWeight','bold','Color','k')
zlabel('Altitude (km)','FontSize',16,'FontWeight','bold','Color','k')

xticks(x); set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
yticks(y); set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
zticks(z); set(gca,'zticklabel',num2str(get(gca,'ztick')','%.2f'))
axis tight

SetFDefaults();
hold off
end