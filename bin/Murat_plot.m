function Murat                  =   Murat_plot(Murat)
%% MURAT_PLOT Creates files for visualization in Matlab and Paraview
%% PATHS and FIGURES
% Importing all the necessary inputs and data for plotting
FLabel                          =   Murat.input.label;
origin1                         =   Murat.input.origin;
ending1                         =   Murat.input.end;
x                               =   Murat.input.x;
y                               =   Murat.input.y;
z                               =   Murat.input.z;
visib                           =   Murat.input.visibility;
sz                              =   Murat.input.sizeMarker;
fformat                         =   Murat.input.format;
sections                        =   Murat.input.sections;
plotV                           =   Murat.input.modvPlot;
Qm                              =   Murat.data.inverseQc;
time0                           =   Murat.data.travelTime;
retainPeakDelay                 =   Murat.data.retainPeakDelay;
retainQc                        =   Murat.data.retainQc;
retainQ                         =   Murat.data.retainQ;
fitrobust                       =   Murat.data.fitrobust;
peakData                        =   Murat.data.peakDelay;
luntot                          =   Murat.data.totalLengthRay;
rma                             =   Murat.data.plotRays;
modv_pd                         =   Murat.data.modvPeakDelay;
modv_Qc                         =   Murat.data.modvQc;
modv_Q                          =   Murat.data.modvQ;
evestazMeters                   =   Murat.data.locationsM;
evestazDegrees                  =   Murat.data.locationsDeg;
sizeTitle                       =   18;
%% 
% Due to the input (lat/long), the code needs to switch 1 and 2. The same happens 
% with events and stations.
FPath                           =   './';
origin                          =   [origin1(2) origin1(1) origin1(3)];
ending                          =   [ending1(2) ending1(1) ending1(3)];
evestaz                         =...
    [evestazDegrees(:,2) evestazDegrees(:,1) -evestazDegrees(:,3)/1000 ...
    evestazDegrees(:,5) evestazDegrees(:,4) evestazDegrees(:,6)/1000];
%% PLOTS - coverage and sensitivity
% Murat_plot starts plotting the ray distribution if asked by the user. It stores 
% the files in the corresponding folder.
storeFolder                     =   'Rays_Checks';
%% 
% The rays are visualized for different techniques, starting with the peak delay
FName_peakDelay             =   'Rays-PeakDelay';
rma_pd                      =   rma(:,2:4,retainPeakDelay)/1000;
evestaz_pd                  =   evestaz(retainPeakDelay,:);
rays_peakDelay              =...
    Murat_imageRays(rma_pd,origin,ending,evestaz_pd,sz,...
    x,y,z,FName_peakDelay,visib);
saveas(rays_peakDelay,...
    fullfile(FPath, FLabel, storeFolder, FName_peakDelay), fformat);
%% 
% The next figure shows the rays for the total attenuation (Q)
FName_Q                     =   'Rays-Q';
rma_Q                       =   rma(:,2:4,retainQ)/1000;
evestaz_Q                   =   evestaz(retainQ,:);
rays_Q                      =...
    Murat_imageRays(rma_Q,origin,ending,evestaz_Q,sz,...
    x,y,z,FName_Q,visib);
saveas(rays_Q,fullfile(FPath, FLabel, storeFolder, FName_Q), fformat);
%% 
% The next sensitivity to check is the one for coda attenuation. The code creates 
% a figure that shows sections in the sensitivity kernels. The left panel shows 
% the sensitivity kernel in the full space while the second shows the normalized 
% kernel in the inversion grid.
FName_Qc                    =   'Kernel-Qc';
kernels                     =...
    figure('Name','Kernels Coda Attenuation',...
    'NumberTitle','off','visible',visib,'Position',[20,400,1200,1000]);
Murat_codaMatrix(Murat,1,...
    Murat.input.startLapseTime,evestazMeters(1,:));
saveas(kernels,fullfile(FPath, FLabel, storeFolder, FName_Qc));
%% Plot - Checks
% In this section Murat_plot checks that (1) Qc is constant with ray length 
% and (2) peak delays increase with travel time.
% 
% It defines the the measurements that have to be kept, then creates the figure.
luntot_Qc                   =   luntot(retainQc)/1000;
mQm                         =   mean(Qm(retainQc));
log10Time                   =   log10(time0(retainPeakDelay));
Qcpd                        =   figure('Name','Qc and peak-delays',...
    'NumberTitle','off','visible',visib,'Position',[20,400,1200,1000]);
%% 
% It starts with Qc relative to the ray length, in the first subplot.
subplot(2,1,1)
plot(luntot_Qc,Qm(retainQc),'o',...
    'MarkerSize',6,'MarkerEdgeColor',[0 0 0])
hold on
plot(luntot_Qc,mQm*ones(length(luntot_Qc),1),'r-','LineWidth',2)
title('Dependence of Qc^{-1} on ray lengths');
xlabel('Ray length (km)','FontSize',sizeTitle,...
    'FontWeight','bold','Color','k')
ylabel('Qc^{-1}','FontSize',sizeTitle,...
    'FontWeight','bold','Color','k')
legend({'Qc^{-1}',cat(2,'<Qc> = ',num2str(1/mQm))},...
    'Location','northeast')
SetFDefaults()
%% 
% Then it shows the peak delay relative to the travel time, in the second subplot.
subplot(2,1,2)
plot(fitrobust,'r-',log10Time,...
    log10(peakData(retainPeakDelay)),'ko')
xti                         =   xticks;
xtl                         =   length(xti);
xt                          =   cell(xtl,1);
for i=1:length(xti)
    xt(i,1)                 =   {10^xti(i)};
end
%% 
% After removing outliers, its shows the best fit parameters for the distance 
% dependence.
nameText                    =   cat(2,'A(f) = ',...
    num2str(fitrobust.p2),', B(f) = ',num2str(fitrobust.p1));
%% 
% The rest is just figure creation and labelling.
xticklabels(xt)
title('Logarithmic dependence of peak delays on travel times');
xlabel('Travel time (s)','FontSize',sizeTitle,'FontWeight','bold',...
    'Color','k')
ylabel('Log. Peak delay','FontSize',sizeTitle,...
    'FontWeight','bold','Color','k')
legend('Location','southeast')
text(min(log10Time) + (mean(log10Time)-min(log10Time))/3,...
    max(log10(peakData)),...
    nameText,'HorizontalAlignment','center','FontSize',20,'Color','r')
SetFDefaults()
FName                       =   'Qc_Peak_Delay';
saveas(Qcpd, fullfile(FPath, FLabel, storeFolder, FName), fformat);
%% PLOT - RESULTS
% Set up matrices. The points are set to the upper SW vertices so that they 
% work with the function "slice". All stored in the appropriate sub-folder.
storeFolder                     =   'Results';
[X,Y,Z1,mPD]                    =   Murat_fold(x,y,z,modv_pd(:,4));
[~,~,~,mQc]                     =   Murat_fold(x,y,z,modv_Qc(:,4));
[~,~,~,mQ]                      =   Murat_fold(x,y,z,modv_Q(:,4));
Z                               =   Z1/1000;
sections(3)                     =   sections(3)/1000;
%% 
% As usual need to select the retained measurements. In addition, use a specific 
% function to generate different color maps.
evestaz_Qc                      =   evestaz(retainQc,:);
cyanpink                        =   colMapGen([1,0,1],[0,1,1],256);
purpleorange                    =...
    colMapGen([0.5 0 0.5],[0.91 0.41 0.17],256);
%% 
% Peak delays results
FName_PDMap                 =   'Peak-Delay-3D';
peakDelaymap                =...
    Murat_image3D(X,Y,Z,mPD,...
    redblue,sections,evestaz_pd,x,y,z,FName_PDMap,sz,visib);
title('Peak-delay variations',...
    'FontSize',sizeTitle,'FontWeight','bold','Color','k');
saveas(peakDelaymap,fullfile(FPath, FLabel, storeFolder, FName_PDMap));
%% 
% Qc results
FName_QcMap                 =   'Qc-3D';
Qcmap                       =...
    Murat_image3D(X,Y,Z,mQc,...
    cyanpink,sections,evestaz_Qc,x,y,z,FName_QcMap,sz,visib);
title('Coda attenuation',...
    'FontSize',sizeTitle,'FontWeight','bold','Color','k');
saveas(Qcmap,fullfile(FPath, FLabel, storeFolder, FName_QcMap));
%% 
% Q results
FName_QMap                  =   'Q-3D';
Qmap                        =...
    Murat_image3D(X,Y,Z,mQ,...
    purpleorange,sections,evestaz_Q,x,y,z,FName_QMap,sz,visib);
title('Total attenuation',...
    'FontSize',sizeTitle,'FontWeight','bold','Color','k');
saveas(Qmap,fullfile(FPath, FLabel, storeFolder, FName_QMap));
%% 
% Also showing the velocity model in case it is available
if Murat.input.availableVelocity == 1
    FName_Vimage                    =   'Velocity_model';
    Vimage                          =   Murat_image3D(X,Y,Z,plotV,...
        inferno,sections,evestaz_Q,x,y,z,FName_Vimage,sz,visib);
    title('Velocity Model',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    
    saveas(Vimage,fullfile(FPath, FLabel, storeFolder, FName_Vimage));
end
%% PLOT - CHECKERBOARDS
% In this section Murat_plot shows input and output of the checkerboard tests 
% for Q and Qc.
[~,~,~,check_inputQc]       =   Murat_fold(x,y,z,modv_Qc(:,6));
[~,~,~,check_outputQc]      =   Murat_fold(x,y,z,modv_Qc(:,7));
[~,~,~,check_inputQ]        =   Murat_fold(x,y,z,modv_Qc(:,6));
[~,~,~,check_outputQ]       =   Murat_fold(x,y,z,modv_Q(:,7));
%% 
% Checkerboard Qc: Input and Output
FNameQcCheck                =   'Qc-Checkerboard';
Qc_check                    =   figure('Name',FNameQcCheck,...
    'NumberTitle','off','visible',visib,'Position',[20,400,2000,1000]);
subplot(1,2,1)
Murat_image3D_2panels(X,Y,Z,check_inputQc,...
    'bone',sections,evestaz_Qc,x,y,z,sz);
title('Input checkerboard Qc',...
    'FontSize',sizeTitle,'FontWeight','bold','Color','k');
subplot(1,2,2)
Murat_image3D_2panels(X,Y,Z,check_outputQc,...
    'bone',sections,evestaz_Qc,x,y,z,sz);
title('Output checkerboard Qc',...
    'FontSize',sizeTitle,'FontWeight','bold','Color','k');
saveas(Qc_check,fullfile(FPath, FLabel, storeFolder, FNameQcCheck));
%Checkerboard Q: Input and Output
FNameQCheck                 =   'Q-Checkerboard';
Q_check                     =   figure('Name',FNameQCheck,...
    'NumberTitle','off','visible',visib,'Position',[20,400,2000,1000]);
subplot(1,2,1)
Murat_image3D_2panels(X,Y,Z,check_inputQ,...
    'bone',sections,evestaz_Q,x,y,z,sz);
title('Input checkerboard Q',...
    'FontSize',sizeTitle,'FontWeight','bold','Color','k');
subplot(1,2,2)
Murat_image3D_2panels(X,Y,Z,check_outputQ,...
    'bone',sections,evestaz_Q,x,y,z,sz);
title('Output checkerboard Q',...
    'FontSize',sizeTitle,'FontWeight','bold','Color','k');
saveas(Q_check,fullfile(FPath, FLabel, storeFolder, FNameQCheck));
%% PLOT - SPIKES
% In this section Murat_plot shows input and output of the spike tests for
% Q and Qc.
[~,~,~,spike_inputQc]       =   Murat_fold(x,y,z,modv_Qc(:,8));
[~,~,~,spike_outputQc]      =   Murat_fold(x,y,z,modv_Qc(:,9));
[~,~,~,spike_inputQ]        =   Murat_fold(x,y,z,modv_Qc(:,8));
[~,~,~,spike_outputQ]       =   Murat_fold(x,y,z,modv_Q(:,9));
%% 
% Spike Qc: Input and Output
FNameQcSpike                =   'Qc-Spike';
Qc_spike                    =   figure('Name',FNameQcSpike,...
    'NumberTitle','off','visible',visib,'Position',[20,400,2000,1000]);
subplot(1,2,1)
Murat_image3D_2panels(X,Y,Z,spike_inputQc,...
    cyanpink,sections,evestaz_Qc,x,y,z,sz);
title('Input spike Qc',...
    'FontSize',sizeTitle,'FontWeight','bold','Color','k');
subplot(1,2,2)
Murat_image3D_2panels(X,Y,Z,spike_outputQc,...
    cyanpink,sections,evestaz_Qc,x,y,z,sz);
title('Output spike Qc',...
    'FontSize',sizeTitle,'FontWeight','bold','Color','k');
saveas(Qc_spike,fullfile(FPath, FLabel, storeFolder, FNameQcSpike));
%% 
% Spike Q: Input and Output
FNameQSpike                 =   'Q-Spike-Input';
Q_spike                     =   figure('Name',FNameQSpike,...
    'NumberTitle','off','visible',visib,'Position',[20,400,2000,1000]);
subplot(1,2,1)
Murat_image3D_2panels(X,Y,Z,spike_inputQ,...
    purpleorange,sections,evestaz_Q,x,y,z,sz);
title('Input spike Q',...
    'FontSize',sizeTitle,'FontWeight','bold','Color','k');
subplot(1,2,2)
Murat_image3D_2panels(X,Y,Z,spike_outputQ,...
    purpleorange,sections,evestaz_Q,x,y,z,sz);
title('Output spike Q',...
    'FontSize',sizeTitle,'FontWeight','bold','Color','k');
saveas(Q_spike,fullfile(FPath, FLabel, storeFolder, FNameQSpike));
%% PARAMETER PLOT
% The final figure output is the parameter plot separation.
storeFolder                     =   'Results';
%% 
% First Qc and Peak delay are separated in 4 quadrants in their parameter space.
[param_plot,par,para_map]   =...
    Murat_imageParameters(x,y,z,modv_pd,modv_Qc,visib);
xlabel('Qc','FontSize',sizeTitle,'FontWeight','bold','Color','k')
ylabel('Log. peak delay','FontSize',sizeTitle,...
    'FontWeight','bold','Color','k')
title('Parameter space plot',...
    'FontSize',sizeTitle,'FontWeight','bold','Color','k');
axis square
FName_Parameters            =   'Parameter_space_variations';
saveas(param_plot,...
    fullfile(FPath, FLabel, storeFolder, FName_Parameters), fformat);
%% 
% The second part produces the spatial plot, setting each node to the corresponding 
% color.
for k=1:length(par(:,1))
    locate_para             =   para_map(:,1) == par(k,1) &...
        para_map(:,2) == par(k,2) & para_map(:,3) == par(k,3);
    para_map(locate_para,4) =   par(k,4);
end
[~,~,~,para_3D]             =   Murat_fold(x,y,z,para_map(:,4)-5);
un_X = unique(para_3D);
fu                          =   find(un_X==-1,1);
%% 
% The four options are: (1) high for both (red); (2) low for both (green); (3) 
% high for peak delays only (cyan); (4) high for inverse Qc only (orange).
if isempty(fu)
    
    cma_para                =...
        [0.7 0.7 0.7;  0 0.8 0; 0 0.6 1; 1 0.6 0];
    HTick                   =...
        {'Average','Ls La','Hs La','Ls Ha'};
    
else
    
    cma_para                =...
        [0.7 0.7 0.7;  0 0.8 0; 0 0.6 1; 1 0.6 0; 1 0 0];
    HTick                   =...
        {'Average','Ls La','Hs La','Ls Ha','Hs Ha'};
    
end
FName_PMap                  =   'Parameter-Map';
ParaMap                     =...
    Murat_image3D(X,Y,Z,para_3D,...
    'bone',sections,evestaz_Qc,x,y,z,FName_Parameters,sz,visib);
colormap(cma_para)
colorbar('Ticks',un_X,'TickLabels',HTick);
title('Parameter separation map',...
    'FontSize',sizeTitle,'FontWeight','bold','Color','k');
saveas(ParaMap,fullfile(FPath, FLabel, storeFolder, FName_PMap));
%% SAVE all results as VTK for visualization in PARAVIEW
storeFolder                     =   'VTK';
[WE_origin, SN_origin]          =   deg2utm(origin(2),origin(1));
%% 
% Converting Lon/Lat to km for paraview visualization with ndgrid
x_origin                        =   x - origin(1);
y_origin                        =   y - origin(2);
UTM_WE                          =   WE_origin + deg2km(x_origin)*1000;
UTM_SN                          =   SN_origin + deg2km(y_origin)*1000;
[X_UTM,Y_UTM,~]                 =   ndgrid(UTM_WE,UTM_SN,z);
%% 
% Writes the four models to vtk
vtkwrite(fullfile(FPath, FLabel,storeFolder,'Peak_delay.vtk'),...
    'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Peak_delay',mPD)
vtkwrite(fullfile(FPath, FLabel,storeFolder,'Qc.vtk'),'structured_grid',...
    X_UTM,Y_UTM,Z1,'scalars','Qc',mQc)
vtkwrite(fullfile(FPath, FLabel,storeFolder,'Q.vtk'),'structured_grid',...
    X_UTM,Y_UTM,Z1,'scalars','Q',mQ)
vtkwrite(fullfile(FPath, FLabel,storeFolder,'Velocity_model.vtk'),...
    'structured_grid',X_UTM,Y_UTM,Z1,'scalars','V',plotV)
%% 
% Write the input-output checkerboard
vtkwrite(fullfile(FPath, FLabel,storeFolder,'Input_check_Qc.vtk'),...
    'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Input_check_Qc',...
    check_inputQc)
vtkwrite(fullfile(FPath, FLabel,storeFolder,'Input_check_Q.vtk'),...
    'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Input_check_Q',...
    check_inputQ)
vtkwrite(fullfile(FPath, FLabel,storeFolder,'Output_check_Qc.vtk'),...
    'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Output_check_Qc',...
    check_outputQc)
vtkwrite(fullfile(FPath, FLabel,storeFolder,'Output_check_Q.vtk'),...
    'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Output_check_Q',...
    check_outputQ)
%% 
% Writes the input-output spikes
vtkwrite(fullfile(FPath, FLabel,storeFolder,'Input_spike_Qc.vtk'),...
    'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Input_spike_Qc',...
    spike_inputQc)
vtkwrite(fullfile(FPath, FLabel,storeFolder,'Input_spike_Q.vtk'),...
    'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Input_spike_Q',...
    spike_inputQ)
vtkwrite(fullfile(FPath, FLabel,storeFolder,'Output_spike_Qc.vtk'),...
    'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Output_spike_Qc',...
    spike_outputQc)
vtkwrite(fullfile(FPath, FLabel,storeFolder,'Output_spike_Q.vtk'),...
    'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Output_spike_Q',...
    spike_outputQ)