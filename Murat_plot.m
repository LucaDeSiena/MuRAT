%%  Creates maps
function Murat                  =   Murat_plot(Murat)


%PATHS and FIGURES
FPath                           =   Murat.input.workingDirectory;
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
figuresData                     =   Murat.input.data;
figuresGeometry                 =   Murat.input.geometry;
figuresInversion                =   Murat.input.inversion;
figuresCheckerboard             =   Murat.input.checkerboard;
figuresSpike                    =   Murat.input.spike;

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

% For plotting needs to switch 1 and 2
origin                          =   [origin1(2) origin1(1) origin1(3)];
ending                          =   [ending1(2) ending1(1) ending1(3)];

% Also evestaz is switched
evestaz                         =...
    [evestazDegrees(:,2) evestazDegrees(:,1) -evestazDegrees(:,3)/1000 ...
    evestazDegrees(:,5) evestazDegrees(:,4) evestazDegrees(:,6)/1000];

% Size of title
sizeTitle                       =   18;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLOTS - coverage and sensitivity

if figuresGeometry == 1
    % Rays for different techniques
    
    % Plotting rays for peak delay
    FName_peakDelay             =   'Rays-PeakDelay';
    
    %Creates figure with rays
    rma_pd                      =   rma(:,2:4,retainPeakDelay)/1000;
    evestaz_pd                  =   evestaz(retainPeakDelay,:);
    
    rays_peakDelay              =...
        Murat_imageRays(rma_pd,origin,ending,evestaz_pd,sz,...
        x,y,z,FName_peakDelay,visib);
    
    saveas(rays_peakDelay,...
        fullfile(FPath, FLabel, FName_peakDelay), fformat);
    
    % Plotting rays for Q
    FName_Q                     =   'Rays-Q';
    
    %Creates figure with rays
    rma_Q                       =   rma(:,2:4,retainQ)/1000;
    evestaz_Q                   =   evestaz(retainQ,:);
    
    rays_Q                      =...
        Murat_imageRays(rma_Q,origin,ending,evestaz_Q,sz,...
        x,y,z,FName_Q,visib);
    
    saveas(rays_Q,fullfile(FPath, FLabel, FName_Q), fformat);
    
    % Plotting kernels for Qc
    FName_Qc                    =   'Kernel-Qc';
    
    kernels                     =...
        figure('Name','Kernels Coda Attenuation',...
        'NumberTitle','off','visible',visib,'Position',[20,400,1200,1000]);
    
    % Showing sensitivity kernels for coda attenuation
    Murat_codaMatrix(Murat,1,...
        Murat.input.startLapseTime,evestazMeters(1,:));
    
    saveas(kernels,fullfile(FPath, FLabel, FName_Qc));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot - Checks
% To check that Qc is constant with ray length and peak delays
% increase with travel time
if figuresData == 1
    luntot_Qc                   =   luntot(retainQc)/1000;
    mQm                         =   mean(Qm(retainQc));
    
    Qcpd                        =   figure('Name','Qc and peak-delays',...
        'NumberTitle','off','visible',visib,'Position',[20,400,1200,1000]);
    
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
    
    log10Time                   =   log10(time0(retainPeakDelay));
    
    subplot(2,1,2)
    plot(fitrobust,'r-',log10Time,...
        log10(peakData(retainPeakDelay)),'ko')
    xti                         =   xticks;
    xtl                         =   length(xti);
    xt                          =   cell(xtl,1);
    for i=1:length(xti)
        xt(i,1)                 =   {10^xti(i)};
    end
    
    nameText                    =   cat(2,'A(f) = ',...
        num2str(fitrobust.p2),', B(f) = ',num2str(fitrobust.p1));
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
    saveas(Qcpd, fullfile(FPath, FLabel, FName), fformat);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLOT - RESULTS
%Set up matrices
%The points are the upper SW vertces so that they work with the function
%slice
[X,Y,Z1,mPD]                    =   Murat_fold(x,y,z,modv_pd(:,4));
[~,~,~,mQc]                     =   Murat_fold(x,y,z,modv_Qc(:,4));
[~,~,~,mQ]                      =   Murat_fold(x,y,z,modv_Q(:,4));
Z                               =   Z1/1000;
sections(3)                     =   sections(3)/1000;
evestaz_Qc                      =   evestaz(retainQc,:);
cyanpink                        =   colMapGen([1,0,1],[0,1,1],256);
purpleorange                    =...
    colMapGen([0.5 0 0.5],[0.91 0.41 0.17],256);

if figuresInversion == 1
    
    
    %Peak delays
    FName_PDMap                 =   'Peak-Delay-3D';
    
    peakDelaymap                =...
        Murat_image3D(X,Y,Z,mPD,...
        redblue,sections,evestaz_pd,x,y,z,FName_PDMap,sz,visib);
    
    title('Peak-delay variations',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    
    saveas(peakDelaymap,fullfile(FPath, FLabel, FName_PDMap));
    
    %Qc
    
    FName_QcMap                 =   'Qc-3D';
    Qcmap                       =...
        Murat_image3D(X,Y,Z,mQc,...
        cyanpink,sections,evestaz_Qc,x,y,z,FName_QcMap,sz,visib);
    
    title('Coda attenuation',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    
    saveas(Qcmap,fullfile(FPath, FLabel, FName_QcMap));
    
    %Q
    purpleorange = colMapGen([.5 0 .5],[0.91 0.41 0.17],256);
    
    FName_QMap                  =   'Q-3D';
    Qmap                        =...
        Murat_image3D(X,Y,Z,mQ,...
        purpleorange,sections,evestaz_Q,x,y,z,FName_QMap,sz,visib);
    
    title('Coda attenuation',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    
    saveas(Qmap,fullfile(FPath, FLabel, FName_QMap));
end
%% PLOT - CHECKERBOARDS

if figuresCheckerboard == 1
    
    [~,~,~,check_inputQc]       =   Murat_fold(x,y,z,modv_Qc(:,6));
    [~,~,~,check_outputQc]      =   Murat_fold(x,y,z,modv_Qc(:,7));
    [~,~,~,check_inputQ]        =   Murat_fold(x,y,z,modv_Qc(:,6));
    [~,~,~,check_outputQ]       =   Murat_fold(x,y,z,modv_Q(:,7));
    
    %Checkerboard Qc
    %Input
    FName_QcInput               =   'Qc-Checkerboard-Input';
    Qcinput                     =...
        Murat_image3D(X,Y,Z,check_inputQc,...
        'bone',sections,evestaz_Qc,x,y,z,FName_QcInput,sz,visib);
    
    title('Input checkerboard Qc',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    
    saveas(Qcinput,fullfile(FPath, FLabel, FName_QcInput));
    
    %Output
    FName_QcOutput              =   'Qc-Checkerboard-Output';
    Qcoutput                    =...
        Murat_image3D(X,Y,Z,check_outputQc,...
        'bone',sections,evestaz_Qc,x,y,z,FName_QcOutput,sz,visib);
    
    title('Output checkerboard Qc',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    
    saveas(Qcoutput,fullfile(FPath, FLabel, FName_QcOutput));
    
    %Checkerboard Q
    %Input
    FName_QInput                =   'Q-Checkerboard-Input';
    Qinput                      =...
        Murat_image3D(X,Y,Z,check_inputQ,...
        'bone',sections,evestaz_Q,x,y,z,FName_QInput,sz,visib);
    
    title('Input checkerboard Q',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    
    saveas(Qinput,fullfile(FPath, FLabel, FName_QInput));
    
    %Output
    FName_QOutput               =   'Q-Checkerboard-Output';
    Qoutput                     =...
        Murat_image3D(X,Y,Z,check_outputQ,...
        'bone',sections,evestaz_Q,x,y,z,FName_QOutput,sz,visib);
    
    title('Output checkerboard Q',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    
    saveas(Qoutput,fullfile(FPath, FLabel, FName_QOutput));
end

%% PLOT - SPIKES

if figuresSpike == 1
    
    [~,~,~,spike_inputQc]       =   Murat_fold(x,y,z,modv_Qc(:,8));
    [~,~,~,spike_outputQc]      =   Murat_fold(x,y,z,modv_Qc(:,9));
    [~,~,~,spike_inputQ]        =   Murat_fold(x,y,z,modv_Qc(:,8));
    [~,~,~,spike_outputQ]       =   Murat_fold(x,y,z,modv_Q(:,9));
    
    %Spike Qc
    %Input
    FName_QcInput               =   'Qc-Spike-Input';
    Qcinput                     =...
        Murat_image3D(X,Y,Z,spike_inputQc,...
        cyanpink,sections,evestaz_Qc,x,y,z,FName_QcInput,sz,visib);
    
    title('Input spike Qc',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    
    saveas(Qcinput,fullfile(FPath, FLabel, FName_QcInput));
    
    %Output
    FName_QcOutput              =   'Qc-Spike-Output';
    Qcoutput                    =...
        Murat_image3D(X,Y,Z,spike_outputQc,...
        cyanpink,sections,evestaz_Qc,x,y,z,FName_QcOutput,sz,visib);
    
    title('Output spike Qc',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    
    saveas(Qcoutput,fullfile(FPath, FLabel, FName_QcOutput));
    
    %Spike Q
    %Input
    FName_QInput                =   'Q-Spike-Input';
    Qinput                      =...
        Murat_image3D(X,Y,Z,spike_inputQ,...
        purpleorange,sections,evestaz_Q,x,y,z,FName_QInput,sz,visib);
    
    title('Input spike Q',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    
    saveas(Qinput,fullfile(FPath, FLabel, FName_QInput));
    
    %Output
    FName_QOutput               =   'Q-Spike-Output';
    Qoutput                     =...
        Murat_image3D(X,Y,Z,spike_outputQ,...
        purpleorange,sections,evestaz_Q,x,y,z,FName_QOutput,sz,visib);
    
    title('Output spike Q',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    
    saveas(Qoutput,fullfile(FPath, FLabel, FName_QOutput));
end

%% PARAMETER PLOT
if figuresInversion == 1
    
    [param_plot,par,para_map]   =...
        Murat_imageParameters(x,y,z,modv_pd,modv_Qc,visib);
    xlabel('Qc','FontSize',sizeTitle,'FontWeight','bold','Color','k')
    ylabel('Log. peak delay','FontSize',sizeTitle,...
        'FontWeight','bold','Color','k')
    title('Parameter space plot',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    axis square
    
    FName_Parameters            =   'Parameter_space_variations';
    saveas(param_plot,fullfile(FPath, FLabel, FName_Parameters), fformat);
    
    for k=1:length(par(:,1))
        locate_para             =   para_map(:,1) == par(k,1) &...
            para_map(:,2) == par(k,2) & para_map(:,3) == par(k,3);
        para_map(locate_para,4) =   par(k,4);
    end
    
    [~,~,~,para_3D]             =   Murat_fold(x,y,z,para_map(:,4)-5);
    un_X = unique(para_3D);
    fu                          =   find(un_X==-1,1);
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
        'bone',sections,evestaz_Qc,x,y,z,FName_QcInput,sz,visib);
    
    colormap(cma_para)
    colorbar('Ticks',un_X,'TickLabels',HTick);
    title('Parameter separation map',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    
    saveas(ParaMap,fullfile(FPath, FLabel, FName_PMap));
end

%% SAVE as VTK for visualization in PARAVIEW
[WE_origin, SN_origin]          =   deg2utm(origin(2),origin(1));
% convert Lon/Lat to UTM
x_origin                        =   x - origin(1);
y_origin                        =   y - origin(2);
UTM_WE                          =   WE_origin + deg2km(x_origin)*1000;
UTM_SN                          =   SN_origin + deg2km(y_origin)*1000;
[X_UTM,Y_UTM,~]                 =   Murat_fold(UTM_WE,UTM_SN,z);

%write the models to vtk
vtkwrite(fullfile(FPath, FLabel,'Peak_delay.vtk'),'structured_grid',...
    X_UTM,Y_UTM,Z1,'scalars','Peak_delay',mPD)

vtkwrite(fullfile(FPath, FLabel,'Qc.vtk'),'structured_grid',...
    X_UTM,Y_UTM,Z1,'scalars','Qc',mQc)

vtkwrite(fullfile(FPath, FLabel,'Q.vtk'),'structured_grid',...
    X_UTM,Y_UTM,Z1,'scalars','Q',mQ)

%write the input-output checks
vtkwrite(fullfile(FPath, FLabel,'Input_check_Qc.vtk'),'structured_grid',...
    X_UTM,Y_UTM,Z1,'scalars','Input_check_Qc',check_inputQc)

vtkwrite(fullfile(FPath, FLabel,'Input_check_Q.vtk'),'structured_grid',...
    X_UTM,Y_UTM,Z1,'scalars','Input_check_Q',check_inputQ)

vtkwrite(fullfile(FPath, FLabel,'Output_check_Qc.vtk'),...
    'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Output_check_Qc',...
    check_outputQc)

vtkwrite(fullfile(FPath, FLabel,'Output_check_Q.vtk'),'structured_grid',...
    X_UTM,Y_UTM,Z1,'scalars','Output_check_Q',check_outputQ)

