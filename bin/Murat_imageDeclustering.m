function clustering       =   Murat_imageDeclustering(locDegOriginal,...
    locDegrees,origin,ending,name)
% function clustering       =...
%     Murat_imageDeclustering(locDegOriginal,locDegrees,origin,ending)
%
% PLOTS 2D rays before and after declustering.
%
% Input parameters:
%    locDegOriginal:        original locations
%    locDegrees:            locations after declustering
%    origin:                start of the grid
%    ending:                end of the grid
%    name:                  title of the figure
%
% Output parameters:
%    clustering:            image produced

clustering              =   figure('Name',name,'NumberTitle','off',...
    'Position',[20,400,1200,1000],'visible','off');

plot([locDegOriginal(:,2),locDegOriginal(:,5)],...
        [locDegOriginal(:,1),locDegOriginal(:,4)],'-k','LineWidth',2)
hold on

plot([locDegrees(:,2), locDegrees(:,5)],...
        [locDegrees(:,1),locDegrees(:,4)],'-r','LineWidth',2)
    
scatter(locDegOriginal(:,2),locDegOriginal(:,1),80,'c','MarkerEdgeColor','b',...
    'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

scatter(locDegrees(:,5),locDegrees(:,4),80,'^','MarkerEdgeColor','m',...
    'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)

legend('Original rays','','Rays after declustering','','Earthquakes','Stations')
title('Comparison between original and declustered data')
axis equal
xlim([origin(2) ending(2)]);
ylim([origin(1) ending(1)]);
grid on
xlabel('Longitude (°)','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Latitude (°)','FontSize',16,'FontWeight','bold','Color','k')
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
xtickangle(45)
ytickangle(45)
SetFDefaults();
hold off

end




