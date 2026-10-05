function Murat_hitmap(CoordGrid,rayCrossing,sections,evestaz,name)
%
% HITMAP produces horizontal and vertical cross sections displaying the
%                   seismic ray count as hitmap.
%
% Input parameters:
%    CoordGrid:     3D coordinates
%    rayCrossing:   vector containing number of rays per cell
%

x = round(CoordGrid(:,1),2);
y = round(CoordGrid(:,2),2);
z = round(CoordGrid(:,3),2);
rayCrossing = rayCrossing(:);

A   = [x y -z];
A2 = sortrows(A,{'ascend','ascend','ascend'});
x   = A2(:,1);
y   = A2(:,2);
z   = -A2(:,3);

xu = unique(x);
yu = unique(y);
zu = sort(unique(z),'descend');

nx = numel(xu);
ny = numel(yu);
nz = numel(zu);


[Xp,Yp,Zp]  =   meshgrid(xu,yu,zu);

R = permute(reshape(rayCrossing, [nz, ny, nx]),[2 3 1]);

% ---------------- Horizontal section ----------------
myfig();
imh = slice(Xp, Yp, Zp, R, [], [], -8000);
colormap(hot)
shading flat
set(imh,'EdgeColor','none')
hcb = colorbar;
hcb.FontSize = 18;
caxis([min(rayCrossing) 100])

hold on
scatter3(evestaz(:,2),evestaz(:,1),evestaz(:,3),30,'c',...
    'MarkerEdgeColor','b', 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

scatter3(evestaz(:,5),evestaz(:,4),evestaz(:,6),30,'^',...
    'MarkerEdgeColor','m', 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

xlabel('Longitude (°)','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Latitude (°)','FontSize',16,'FontWeight','bold','Color','k')

xticks(xu); set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
yticks(yu); set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))

view(2)
axis tight
SetFDefaults();
saveFigureAsImage([name '_H']);
close

% ---------------- W-E section ----------------
myfig();
imwe = slice(Xp, Yp, Zp, R, sections(2), [], []);
colormap(hot)
shading flat
set(imwe,'EdgeColor','none')
hcb = colorbar;
hcb.FontSize = 18;
caxis([min(rayCrossing) 100])

hold on
scatter3(evestaz(:,2),evestaz(:,1),evestaz(:,3),30,'c',...
    'MarkerEdgeColor','b', 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

scatter3(evestaz(:,5),evestaz(:,4),evestaz(:,6),30,'^',...
    'MarkerEdgeColor','m', 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

xlabel('Longitude (°)','FontSize',16,'FontWeight','bold','Color','k')
zlabel('Altitude (km)','FontSize',16,'FontWeight','bold','Color','k')

xticks(xu); set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
view(90,0)

axis tight
SetFDefaults();
saveFigureAsImage([name '_WE']);
close

% ---------------- S-N section ----------------
myfig();
imsn = slice(Xp, Yp, Zp, R, [], sections(1), []);
colormap(hot)
shading flat
set(imsn,'EdgeColor','none')
hcb = colorbar;
hcb.FontSize = 18;
caxis([min(rayCrossing) 100])

hold on
scatter3(evestaz(:,2),evestaz(:,1),evestaz(:,3),30,'c',...
    'MarkerEdgeColor','b', 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

scatter3(evestaz(:,5),evestaz(:,4),evestaz(:,6),30,'^',...
    'MarkerEdgeColor','m', 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

ylabel('Latitude (°)','FontSize',16,'FontWeight','bold','Color','k')
zlabel('Altitude (km)','FontSize',16,'FontWeight','bold','Color','k')

yticks(yu); set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
view(0,0)

axis tight
SetFDefaults();
saveFigureAsImage([name '_SN']);
close

end
