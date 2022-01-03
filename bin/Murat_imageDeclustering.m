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
for n = 1:length(locDegOriginal(:,1))
    plot([locDegOriginal(n,2),locDegOriginal(n,5)],...
        [locDegOriginal(n,1),locDegOriginal(n,4)],'-k')
end
for n = 1:length(locDegrees(:,1))
    plot([locDegrees(n,2), locDegrees(n,5)],...
        [locDegrees(n,1),locDegrees(n,4)],'-r')
end

hold on
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




