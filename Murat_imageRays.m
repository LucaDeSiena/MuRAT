function rays           =...
    Murat_imageRays(rma,origin,ending,evestaz,sz,x,y,z,name,visib)
%PLOTS 3D rays

rays                    =   figure('Name',name,...
    'NumberTitle','off','visible',visib,'Position',[20,400,1200,1000]);

load coastlines coastlat coastlon
z                       =   sort(z)/1000;
lrma                    =   length(rma(1,1,:));

for i = 1:lrma
    ray                 =   rma(:,:,i);
    ray(:,1)            =   origin(1) + km2deg(ray(:,1));
    ray(:,2)            =   origin(2) + km2deg(ray(:,2));
    
    subplot(2,2,1)
    hold on
    plot(ray(:,1),ray(:,2),'k');
    geoshow(coastlat,coastlon);
    hold off
    
    subplot(2,2,2)
    hold on
    plot(ray(:,2),ray(:,3),'k');
    hold off
    
    subplot(2,2,3)
    hold on
    plot(ray(:,1),ray(:,3),'k');
    hold off
end


subplot(2,2,1)
hold on
scatter(evestaz(:,1),evestaz(:,2),sz,...
    'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

scatter(evestaz(:,4),evestaz(:,5),sz,...
    '^','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
xlim([origin(1) ending(1)]);
ylim([origin(2) ending(2)]);

xlabel('WE','FontSize',16,'FontWeight','bold','Color','k')
ylabel('SN','FontSize',16,'FontWeight','bold','Color','k')

xticks(x)
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
yticks(y)
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))

SetFDefaults();
hold off

subplot(2,2,2)
hold on
scatter(evestaz(:,2),evestaz(:,3),sz,...
    'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

scatter(evestaz(:,5),evestaz(:,6),sz,...
    '^','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
xlim([origin(2) ending(2)]);

xlabel('SN','FontSize',12,...
    'FontWeight','bold','Color','k')
ylabel('Depth (km)','FontSize',12,'FontWeight','bold','Color','k')
xticks(y)
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
yticks(z)
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
SetFDefaults();
hold off

subplot(2,2,3)
hold on
scatter(evestaz(:,1),evestaz(:,3),sz,...
    'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

scatter(evestaz(:,4),evestaz(:,6),sz,...
    '^','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
xlim([origin(1) ending(1)]);
xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Depth (km)','FontSize',12,...
    'FontWeight','bold','Color','k')
xticks(x)
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
yticks(z)
set(gca,'zticklabel',num2str(get(gca,'ytick')','%.2f'))
SetFDefaults();
hold off

subplot(2,2,4)
hold on

geoshow(coastlat,coastlon,'Color','k');
scatter(evestaz(1,4),evestaz(1,5),100,...
    'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

xticks(x(1))
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
yticks(y(1))
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
axis auto
SetFDefaults();

hold off



