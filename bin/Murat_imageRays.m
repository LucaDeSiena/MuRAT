function rays           =...
    Murat_imageRays(rma,origin,ending,evestaz,x,y,z,name)
% function rays           =...
%     Murat_imageRays(rma,origin,ending,evestaz,x,y,z,name)
%
% PLOTS 3D rays
%
% Input parameters:
%    rma:               ray in the format output from ray tracing
%    origin:            origin of the grid
%    ending:            end of the grid
%    evestaz:           locations of events and stations in meters
%    x:                 x vector
%    y:                 y vector
%    z:                 z vector
%    name:              title of the figure
%
% Output parameters:
%    rays:              image produced

rays                    =   figure('Name',name,...
    'NumberTitle','off','Position',[20,400,1200,1000],'visible','off');

load coastlines coastlat coastlon
z1                      =   sort(z)/1000;
lrma                    =   length(rma(1,1,:));
evestaz_ray             =   zeros(lrma,6);

wgs84                               =   wgs84Ellipsoid("m");

subplot(2,2,1)

for i = 1:lrma
    ray                     =   rma(:,:,i)*1000;
    ray(:,3)                =   rma(:,3,i);
    
    d                       =   sqrt(ray(:,1).^2 + ray(:,2).^2);
    az                      =   atan(ray(:,1)./ray(:,2))*360/2/pi;
    az(isnan(az))           =   0;
    [lat2,lon2]             =   reckon(origin(1),origin(2),d,az,wgs84);

    subplot(2,2,1)
    hold on
    plot(lon2,lat2,'k');
    
    subplot(2,2,2)
    hold on
    plot(lat2,ray(:,3),'k');
    
    subplot(2,2,3)
    hold on
    plot(lon2,ray(:,3),'k');
    
    % save start and end point of ray
    evestaz_ray(i,:)        =   [lat2(1),lon2(1),ray(1,3),...
        lat2(end),lon2(end),ray(end,3)];
end

centreGrid              =   [origin(2) + (ending(2) - origin(2))/2 ...
    origin(1) + (ending(1) - origin(1))/2];

subplot(2,2,1)
hold on
scatter(evestaz_ray(:,2),evestaz_ray(:,1),60,'c','MarkerEdgeColor','b',...
    'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
scatter(evestaz(:,5),evestaz(:,4),60,'^','MarkerEdgeColor','m',...
    'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
axis equal
xlim([origin(2) ending(2)]);
ylim([origin(1) ending(1)]);
xlabel('Longitude (°)','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Latitude (°)','FontSize',16,'FontWeight','bold','Color','k')
xticks(x(1:2:end-1))
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
yticks(y(1:2:end-1))
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
xtickangle(45)
ytickangle(45)
SetFDefaults();
hold off


subplot(2,2,2)
hold on
scatter(evestaz_ray(:,1),evestaz_ray(:,3),60,'c','MarkerEdgeColor','b',...
    'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
scatter(evestaz(:,4),evestaz(:,6),60,'^','MarkerEdgeColor','m',...
    'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
xlim([origin(1) ending(1)]);

xlabel('Latitude (°)','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Altitude (km)','FontSize',16,'FontWeight','bold','Color','k')
xticks(y(1:2:end-1))
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
yticks(z1(1:end-1))
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
xtickangle(45)
ytickangle(45)
SetFDefaults();
hold off

subplot(2,2,3)
hold on
scatter(evestaz_ray(:,2),evestaz_ray(:,3),60,'c','MarkerEdgeColor','b',...
    'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
scatter(evestaz(:,5),evestaz(:,6),60,'^','MarkerEdgeColor','m',...
    'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
xlim([origin(2) ending(2)]);
xlabel('Longitude (°)','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Altitude (km)','FontSize',16,'FontWeight','bold','Color','k')
xticks(x(1:2:end-1))
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
yticks(z1(1:end-1))
set(gca,'zticklabel',num2str(get(gca,'ytick')','%.2f'))
xtickangle(45)
ytickangle(45)
SetFDefaults();
hold off

subplot(2,2,4)
hold on

geoshow(coastlat,coastlon,'Color','k');
scatter(centreGrid(1),centreGrid(2),100,'c','MarkerEdgeColor','g',...
    'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
xlabel('Longitude (°)','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Latitude (°)','FontSize',16,'FontWeight','bold','Color','k')
xticks(centreGrid(1))
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
yticks(centreGrid(2))
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
axis auto
xtickangle(45)
ytickangle(45)
SetFDefaults();

hold off
end




