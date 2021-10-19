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
    'NumberTitle','off','Position',[20,400,1200,1000]);

load coastlines coastlat coastlon
z                       =   sort(z)/1000;
lrma                    =   length(rma(1,1,:));

subplot(2,2,1)
geoshow(coastlat,coastlon);
rma(:,1,:)              =   origin(2) + km2deg(rma(:,1,:));
rma(:,2,:)              =   origin(1) + km2deg(rma(:,2,:));
for i = 1:lrma
    ray                 =   rma(:,:,i);
    
    subplot(2,2,1)
    hold on
    plot(ray(:,1),ray(:,2),'k');
    
    subplot(2,2,2)
    hold on
    plot(ray(:,2),ray(:,3),'k');
    
    subplot(2,2,3)
    hold on
    plot(ray(:,1),ray(:,3),'k');
end

centreGrid              =   [origin(2) + (ending(2) - origin(2))/2 ...
    origin(1) + (ending(1) - origin(1))/2];

subplot(2,2,1)
hold on
scatter(evestaz(:,2),evestaz(:,1),60,'c','MarkerEdgeColor',[1 1 1],...
    'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
scatter(evestaz(:,5),evestaz(:,4),60,'^','MarkerEdgeColor',[1 1 1],...
    'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
xlim([origin(2) ending(2)]);
ylim([origin(1) ending(1)]);

xlabel('WE','FontSize',16,'FontWeight','bold','Color','k')
ylabel('SN','FontSize',16,'FontWeight','bold','Color','k')

xticks(x(1:2:end-1))
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
yticks(y(1:2:end-1))
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))

SetFDefaults();
hold off

subplot(2,2,2)
hold on
scatter(evestaz(:,1),evestaz(:,3),60,'c','MarkerEdgeColor',[1 1 1],...
    'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
scatter(evestaz(:,4),evestaz(:,6),60,'^','MarkerEdgeColor',[1 1 1],...
    'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
xlim([origin(1) ending(1)]);

xlabel('SN','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Depth (km)','FontSize',12,'FontWeight','bold','Color','k')
xticks(y(1:2:end-1))
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
yticks(z(1:2:end-1))
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
SetFDefaults();
hold off

subplot(2,2,3)
hold on
scatter(evestaz(:,2),evestaz(:,3),60,'c','MarkerEdgeColor',[1 1 1],...
    'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
scatter(evestaz(:,5),evestaz(:,6),60,'^','MarkerEdgeColor',[1 1 1],...
    'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
xlim([origin(2) ending(2)]);
xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Depth (km)','FontSize',12,'FontWeight','bold','Color','k')
xticks(x(1:2:end-1))
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
yticks(z(1:2:end-1))
set(gca,'zticklabel',num2str(get(gca,'ytick')','%.2f'))
SetFDefaults();
hold off

subplot(2,2,4)
hold on

geoshow(coastlat,coastlon,'Color','k');
scatter(centreGrid(1),centreGrid(2),100,'c','MarkerEdgeColor',[1 1 1],...
    'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

xticks(centreGrid(1))
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
yticks(centreGrid(2))
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
axis auto
SetFDefaults();

hold off
end




