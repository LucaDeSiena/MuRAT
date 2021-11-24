function pd_analysis        =   Murat_imageCheckPeakDelay(time0PD,...
    fitrobust_k,peakData_k,sizeTitle,pd_title)
% function pd_analysis        =   Murat_imageCheckPeakDelay(rtpdk,...
%     time0,fitrobust,peakData_k,sizeTitle,pd_title)
%
% PLOTS the Peak Delay check
%
% Input parameters:
%    time0PD:       travel time for peak delay
%    fitrobust:     robustly fitted coefficients from peak delay fit
%    peakData_k:    peak delay data
%    sizeTitle:     size of title font
%    pd_title:      title of peak delay figure
%
% Output parameters:
%    pd_analysis:    figure for peak delay check

pd_analysis                 =   figure('Name',pd_title,...
    'NumberTitle','off','Position',[20,400,1200,1000],'visible','off');

log10Time                   =   log10(time0PD);
fitrobust_i                 =...
    fitrobust_k(1)*log10Time + fitrobust_k(2);

plot(log10Time,fitrobust_i,'r-',log10Time,log10(peakData_k),'ko')

xti                         =   xticks;
xt                          =   cell(length(xti),1);
for i = 1:length(xti)
    xt(i,1)                 =   {10^xti(i)};
end
yti                         =   yticks;
yt                          =   cell(length(yti),1);
for i = 1:length(yticks)
    yt(i,1)                 =   {10^yti(i)};
end
%%
% After removing outliers, it shows the best fit parameters for the
% distance dependence.
xticklabels(xt)
yticklabels(yt)
title(['Dependence of peak delays on travel times: ',...
    'A(f) = ', num2str(fitrobust_k(2)),...
    ', B(f) = ',num2str(fitrobust_k(1))]);
xlabel('Travel time (s)','FontSize',sizeTitle,'FontWeight','bold',...
    'Color','k')
ylabel('Peak delay (s)','FontSize',sizeTitle,...
    'FontWeight','bold','Color','k')
legend('Location','southeast')
SetFDefaults()