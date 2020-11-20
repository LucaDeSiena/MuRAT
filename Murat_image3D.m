function image          = ...
    Murat_image3D(X,Y,Z,V,color,sections,evestaz,x,y,z,name,sz,visib)
%PLOTS a 3D image of a field on slices

stepgx                  =   x(2) - x(1);
stepgy                  =   y(2) - y(1);
stepgz                  =   z(2) - z(1);
divi                    =   5;
divix                   =   stepgx/divi;
diviy                   =   stepgy/divi;
diviz                   =   stepgz/divi;

xp                      =   x(1)-divix:divix:x(end)+divix;
yp                      =   y(1)-diviy:diviy:y(end)+diviy;
zp                      =   z(1)-diviz:diviz:z(end)+diviz;
zp                      =   zp/1000;

[Xp,Yp,Zp]              =   meshgrid(yp,xp,zp);
mVp                     =   interp3(X,Y,Z,V,Xp,Yp,Zp);

mVp                      =   inpaintn(mVp);

z                       =   sort(z)/1000;
% Creates and outputs figure in 3D

image                   =   figure('Name',name,...
    'NumberTitle','off','visible',visib,'Position',[20,400,1200,1000]);

slice(Xp, Yp, Zp, mVp, sections(1), sections(2), sections(3))
colormap(color);
colorbar
shading flat
hcb                     =   colorbar;
hcb.FontSize            =   14;

hold on
scatter3(evestaz(:,2),evestaz(:,1),evestaz(:,3),sz,...
    'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

scatter3(evestaz(:,5),evestaz(:,4),evestaz(:,6),sz,...
    '^','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

xlabel('SN','FontSize',16,'FontWeight','bold','Color','k')
ylabel('WE','FontSize',16,'FontWeight','bold','Color','k')
zlabel('Depth (km)','FontSize',16,'FontWeight','bold','Color','k')

xticks(y)
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
yticks(x)
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
zticks(z)
set(gca,'zticklabel',num2str(get(gca,'ztick')','%.2f'))
axis tight

SetFDefaults();
hold off
