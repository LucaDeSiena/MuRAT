function Murat                  =   Murat_plot(Murat)
%% MURAT_PLOT Creates files for visualization in Matlab and Paraview
% Importing all the necessary inputs and data for plotting
FLabel                          =   Murat.input.label;
origin1                         =   Murat.input.origin;
ending1                         =   Murat.input.end;
x                               =   Murat.input.x;
y                               =   Murat.input.y;
z                               =   Murat.input.z;
fformat                         =   Murat.input.format;
sections                        =   Murat.input.sections;
plotV                           =   Murat.input.modvPlot;
cf                              =   Murat.input.centralFrequency;
nonlinear                       =   Murat.input.nonLinear;
outputLCurve                    =   Murat.input.lCurve;
if outputLCurve == 0
    lCurveQ                     =   Murat.input.lCurveQ;
elseif outputLCurve == 1
    lCurveQ                     =   [];
end

Qm                              =   Murat.data.inverseQc;
time0                           =   Murat.data.travelTime;
retainPeakDelay                 =   Murat.data.retainPeakDelay;
retainQc                        =   Murat.data.retainQc;
retainQ                         =   Murat.data.retainQ;
ray_crosses_pd                  =   Murat.data.raysPeakDelay;
ray_crosses_Qc                  =   Murat.data.raysQc;
ray_crosses_Q                   =   Murat.data.raysQ;
fitrobust                       =   Murat.data.fitrobust;
peakData                        =   Murat.data.peakDelay;
luntot                          =   Murat.data.totalLengthRay;
rma                             =   Murat.data.raysPlot;
modv_pd                         =   Murat.data.modvPeakDelay;
modv_Qc                         =   Murat.data.modvQc;
modv_Q                          =   Murat.data.modvQ;
evestazMeters                   =   Murat.data.locationsM;
evestazDegrees                  =   Murat.data.locationsDeg;
constQc                         =   Murat.data.const_Qc;
energyRatio                     =   Murat.data.energyRatioBodyCoda;
codaNoiseRatio                  =   Murat.data.energyRatioCodaNoise;
Ac_i                            =   Murat.data.inversionMatrixQc;
RZZ                             =   Murat.data.uncertaintyQc;
A_i                             =   Murat.data.inversionMatrixQ;
residualQc                      =   Murat.data.residualQc;
residualQ                       =   Murat.data.residualQ;

FPath                           =   './';
sizeTitle                       =   18;
lMF                             =   size(ray_crosses_pd);
sections(3)                     =   sections(3)/1000;

cyanpink                        =   colMapGen([1,0,1],[0,1,1],256);
purpleorange                    =...
    colMapGen([0.5 0 0.5],[0.91 0.41 0.17],256);

%%
% Due to the input (lat/long), the code needs to switch 1 and 2.
% The same happens with events and stations.
origin                          =   [origin1(2) origin1(1) origin1(3)];
ending                          =   [ending1(2) ending1(1) ending1(3)];
[WE_origin, SN_origin]          =   deg2utm(origin(2),origin(1));
evestaz                         =...
    [evestazDegrees(:,2) evestazDegrees(:,1) -evestazDegrees(:,3)/1000 ...
    evestazDegrees(:,5) evestazDegrees(:,4) evestazDegrees(:,6)/1000];

%% PLOTS - coverage and sensitivity
% Murat_plot starts plotting the ray distribution if asked by the user.
% It stores  the files in the corresponding folder.

for k = 1:lMF(2)
    storeFolder                     =   'Rays_Checks';
    cfk                         =   cf(k);
    fcName                      =   num2str(cfk);
    if find(fcName == '.')
        fcName(fcName == '.')   =   '_';
    end
    rtpdk                       =   retainPeakDelay(:,k);
    rtQk                        =   retainQ(:,k);
    rtQck                       =   retainQc(:,k);
    rcQk                        =   ray_crosses_Q(:,k);
    rcQck                       =   ray_crosses_Qc(:,k);
    %%
    % The rays are visualized for different techniques, starting with the peak delay
    FName_peakDelay             =   ['Rays_PeakDelay_' fcName '_Hz'];
    rma_pd                      =   rma(:,2:4,rtpdk)/1000;
    evestaz_pd                  =   evestaz(rtpdk,:);
    rays_peakDelay              =...
        Murat_imageRays(rma_pd,origin,ending,evestaz_pd,...
        x,y,z,FName_peakDelay);
    saveas(rays_peakDelay,...
        fullfile(FPath, FLabel, storeFolder, FName_peakDelay), fformat);
    close(rays_peakDelay)
    %%
    % The next figure shows the rays for the total attenuation (Q)
    FName_Q                     =   ['Rays_Q_' fcName '_Hz'];
    
    rma_Q                       =   rma(:,2:4,rtQk)/1000;
    evestaz_Q                   =   evestaz(rtQk,:);
    rays_Q                      =...
        Murat_imageRays(rma_Q,origin,ending,evestaz_Q,...
        x,y,z,FName_Q);
    saveas(rays_Q,fullfile(FPath, FLabel, storeFolder, FName_Q), fformat);
    close(rays_Q)
    %%
    % The next sensitivity to check is the one for coda attenuation. The code creates
    % a figure that shows sections in the sensitivity kernels. The left panel shows
    % the sensitivity kernel in the full space while the second shows the normalized
    % kernel in the inversion grid.
    FName_Qc                    =   ['Kernel_Qc' fcName '_Hz'];
    kernels                     =...
        figure('Name',FName_Qc,...
        'NumberTitle','off','Position',[20,400,1200,1000]);
    Murat_codaMatrix(Murat,1,...
        Murat.input.startLapseTime,evestazMeters(1,:));
    saveas(kernels,fullfile(FPath, FLabel, storeFolder, FName_Qc));
    close(kernels)
    
    %% Plot - Checks
    % In this section Murat_plot makes checks on the three parameters.
    % These plots are always visualised. They check that:
    % (1) Qc is constant with ray length;
    % (2) peak delays increase with travel time;
    % (3) amplitude ratios decay with hypocentral distance.
    % These plots are used to select measurements and understand how well
    % they follow the assumptions.
    Qc_title                    =   ['Qc check ' fcName ' Hz'];
    Qm_k                        =   Qm(rtQck,k);
    Qc_analysis                 =   figure('Name',...
        Qc_title,'NumberTitle','off','Position',[20,400,1200,1000]);
    luntot_Qc                   =   luntot(rtQck)/1000;
    mQm                         =   mean(Qm_k);
    Ac                          =   Ac_i(rtQck,rcQck);
    RZZ_k                       =   RZZ(rtQck,k);
    
    Wc                          =   Murat_weighting(nonlinear,RZZ_k);
    Gc                          =   Wc*Ac;
    dck                         =   Wc*Qm_k;
    [Uc,Sc,~]                   =   svd(Gc);
    %%
    % It starts with Qc relative to the ray length in the first plot.
    subplot(2,1,1)
    plot(luntot_Qc,Qm_k,'o',...
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
    % A Picard plot evaluates the quality of the inversion.
    
    subplot(2,1,2)
    picard(Uc,diag(Sc),dck);
    title(['Picard condition. Residual is: ' num2str(residualQc(k))]);
    SetFDefaults()
    FName                       =   ['Qc_analysis_' fcName '_Hz'];
    saveas(Qc_analysis, fullfile(FPath, FLabel, storeFolder, FName),...
        fformat);
    %%
    % Then it shows the peak delay relative to the travel time.
    log10Time                   =   log10(time0(rtpdk));
    fitrobust_i                 =...
        fitrobust(1,k)*log10Time + fitrobust(2,k);
    pd_title                    =   ['Peak Delay check ' fcName ' Hz'];
    pd_analysis                 =   figure('Name',...
        pd_title,'NumberTitle','off','Position',[20,400,1200,1000]);
    plot(log10Time,fitrobust_i,...
        'r-',log10Time,log10(peakData(rtpdk,k)),'ko')
    xti                         =   xticks;
    xtl                         =   length(xti);
    xt                          =   cell(xtl,1);
    for i = 1:length(xti)
        xt(i,1)                 =   {10^xti(i)};
    end
    
    yti                         =   yticks;
    ytl                         =   length(yti);
    yt                          =   cell(ytl,1);
    for i = 1:length(yti)
        yt(i,1)                 =   {10^yti(i)};
    end
    %%
    % After removing outliers, it shows the best fit parameters for the
    % distance dependence.
    xticklabels(xt)
    yticklabels(yt)
    title(['Dependence of peak delays on travel times: ',...
        'A(f) = ', num2str(fitrobust(2,k)),...
        ', B(f) = ',num2str(fitrobust(1,k))]);
    xlabel('Travel time (s)','FontSize',sizeTitle,'FontWeight','bold',...
        'Color','k')
    ylabel('Peak delay (s)','FontSize',sizeTitle,...
        'FontWeight','bold','Color','k')
    legend('Location','southeast')
    SetFDefaults()
    FName                       =   ['PD_analysis_' fcName '_Hz'];
    saveas(pd_analysis, fullfile(FPath, FLabel, storeFolder, FName), fformat);
    %%
    % Then it plots first the logarithm of the energy ratio versus travel
    % time.
    A                           =   A_i(rtQk,rcQk);
    energyRatio_k               =   energyRatio(rtQk,k);
    const_Qc_k                  =   constQc(rtQk,k);
    [~,~,~,~,d1k,spreadAverageQ]=...
                                Murat_tikhonovQ(cfk,rtQk,outputLCurve,...
                                energyRatio_k,const_Qc_k,...
                                luntot,time0,A,lCurveQ,0);
    close
    
    equationQ                   =   -log(constQc(rtQk,k))...
        -spreadAverageQ(1,1)*log(luntot(rtQk))...
        -2*pi*cfk*time0(rtQk)*spreadAverageQ(2,1);
    dRatio                      =	log(energyRatio_k);
    luntotQ                     =   luntot(rtQk)/1000;
    Edirect                     =...
        energyRatio_k./codaNoiseRatio(rtQk,k);
    [U,S,~]                     =   svd(A);
    
    CN_title                    =...
        ['Coda Normalization check ' fcName ' Hz'];
    CN_analysis                 =   figure('Name',...
        CN_title,'NumberTitle','off','Position',[20,400,1200,1000]);
    %%
    %Plot of the left and right hand sides of the CN equation.
    subplot(3,1,1)
    plot(time0(rtQk),dRatio, 'o',...
        'MarkerSize',6,'MarkerEdgeColor',[0 0 0])
    hold on
    plot(time0(rtQk),equationQ,'r*','MarkerSize',6)
    hold off
    xlabel('Travel time (s)','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('Logarithmic direct-to-coda energy ratio','FontSize',12,...
        'FontWeight','bold','Color','k')
    title(['Dependence of logarithmic energy ratios on travel time. '...
        'Geometrical spreading is: ' num2str(spreadAverageQ(1,1))...
        '+/- ' num2str(spreadAverageQ(1,2))]);
    legend({'Q^{-1}',cat(2,'<Q> = ',num2str(1/spreadAverageQ(2,1)))},...
        'Location','northeast')
    SetFDefaults()
    %%
    %Plot of the direct energy versus distance.
    subplot(3,1,2)
    plot(log(luntotQ),log(Edirect),'o','MarkerSize',6,...
        'MarkerEdgeColor',[0 0 0])
    xlabel('Hypocentral distance (km)','FontSize',12,...
        'FontWeight','bold','Color','k')
    ylabel('Logarithmic direct energy (J/m^2)','FontSize',12,...
        'FontWeight','bold','Color','k')
    title('Dependence of logarithmic direct energy on hypocentral distance.');
    xti                             =   xticks;
    xtl                             =   length(xti);
    xt                              =   cell(xtl,1);
    for i=1:length(xti)
        xt(i,1)                     =   {exp(xti(i))};
    end
    xticklabels(xt)
    SetFDefaults()
    %%
    % Then it checks the inversion with a Picard
    subplot(3,1,3)
    picard(U,diag(S),d1k);
    title(['Picard condition. Residual is: ' num2str(residualQ(k))]);
    SetFDefaults()
    FName                       =   ['CN_analysis_' fcName '_Hz'];
    saveas(CN_analysis, fullfile(FPath, FLabel, storeFolder, FName), fformat);
    
    %% PLOT - RESULTS
    % Set up matrices. The points are set to the upper SW vertices to
    % work with the function "slice". All stored in the sub-folder.
    storeFolder                     =   'Results';
    modv_pd_k                       =   modv_pd(:,:,k);
    modv_Qc_k                       =   modv_Qc(:,:,k);
    modv_Q_k                        =   modv_Q(:,:,k);
    [X,Y,Z1,mPD]                    =   Murat_fold(x,y,z,modv_pd_k(:,5));
    [~,~,~,mQc]                     =   Murat_fold(x,y,z,modv_Qc_k(:,5));
    [~,~,~,mQ]                      =   Murat_fold(x,y,z,modv_Q_k(:,5));
    Z                               =   Z1/1000;
    %%
    % Need to select the retained measurements. Use a specific
    % function to generate different color maps.
    evestaz_Qc                      =   evestaz(rtQck,:);

    %%
    % Peak delays results
    FName_PDMap                     =   ['Peak-Delay-3D_' fcName '_Hz'];
    peakDelaymap                    =...
        Murat_image3D(X,Y,Z,mPD,...
        redblue,sections,evestaz_pd,x,y,z,FName_PDMap);
    title('Peak-delay variations',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    saveas(peakDelaymap,fullfile(FPath, FLabel, storeFolder, FName_PDMap));
    close(peakDelaymap)

    %%
    % Qc results
    FName_QcMap                     =   ['Qc-3D_' fcName '_Hz'];
    Qcmap                           =...
        Murat_image3D(X,Y,Z,mQc,...
        cyanpink,sections,evestaz_Qc,x,y,z,FName_QcMap);
    title('Coda attenuation',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    saveas(Qcmap,fullfile(FPath, FLabel, storeFolder, FName_QcMap));
    close(Qcmap)

    %%
    % Q results
    FName_QMap                      =   ['Q-3D_' fcName '_Hz'];
    Qmap                            =...
        Murat_image3D(X,Y,Z,mQ,...
        purpleorange,sections,evestaz_Q,x,y,z,FName_QMap);
    title('Total attenuation',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    saveas(Qmap,fullfile(FPath, FLabel, storeFolder, FName_QMap));
    close(Qmap)
    
    %% PLOT - CHECKERBOARDS
    % In this section Murat_plot shows the checkerboard tests
    % for Q and Qc.
    storeFolder                     =   'Resolution';
    [~,~,~,check_inputQc]           =   Murat_fold(x,y,z,modv_Qc_k(:,6));
    [~,~,~,check_outputQc]          =   Murat_fold(x,y,z,modv_Qc_k(:,7));
    [~,~,~,check_inputQ]            =   Murat_fold(x,y,z,modv_Qc_k(:,6));
    [~,~,~,check_outputQ]           =   Murat_fold(x,y,z,modv_Q_k(:,7));
    %%
    % Checkerboard Qc: Input and Output
    FNameQcCheck                    =   ['Qc-Checkerboard_' fcName '_Hz'];
    Qc_check                        =   figure('Name',FNameQcCheck,...
        'NumberTitle','off','Position',[20,400,2000,1000]);
    subplot(1,2,1)
    Murat_image3D_2panels(X,Y,Z,check_inputQc,...
        'bone',sections,evestaz_Qc,x,y,z);
    title('Input checkerboard Qc',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    subplot(1,2,2)
    Murat_image3D_2panels(X,Y,Z,check_outputQc,...
        'bone',sections,evestaz_Qc,x,y,z);
    title('Output checkerboard Qc',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    saveas(Qc_check,fullfile(FPath, FLabel, storeFolder, FNameQcCheck));
    close(Qc_check)
    
    %%
    %Checkerboard Q: Input and Output
    FNameQCheck                     =   ['Q-Checkerboard_' fcName '_Hz'];
    Q_check                         =   figure('Name',FNameQCheck,...
        'NumberTitle','off','Position',[20,400,2000,1000]);
    subplot(1,2,1)
    Murat_image3D_2panels(X,Y,Z,check_inputQ,...
        'bone',sections,evestaz_Q,x,y,z);
    title('Input checkerboard Q',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    subplot(1,2,2)
    Murat_image3D_2panels(X,Y,Z,check_outputQ,...
        'bone',sections,evestaz_Q,x,y,z);
    title('Output checkerboard Q',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    saveas(Q_check,fullfile(FPath, FLabel, storeFolder, FNameQCheck));
    close(Q_check)
    
    %% PLOT - SPIKES
    % In this section Murat_plot shows input and output of the spike tests for
    % Q and Qc.
    [~,~,~,spike_inputQc]           =   Murat_fold(x,y,z,modv_Qc_k(:,8));
    [~,~,~,spike_outputQc]          =   Murat_fold(x,y,z,modv_Qc_k(:,9));
    [~,~,~,spike_inputQ]            =   Murat_fold(x,y,z,modv_Qc_k(:,8));
    [~,~,~,spike_outputQ]           =   Murat_fold(x,y,z,modv_Q_k(:,9));
    %%
    % Spike Qc: Input and Output
    FNameQcSpike                    =   ['Qc-Spike_' fcName '_Hz'];
    Qc_spike                        =   figure('Name',FNameQcSpike,...
        'NumberTitle','off','Position',[20,400,2000,1000]);
    subplot(1,2,1)
    Murat_image3D_2panels(X,Y,Z,spike_inputQc,...
        cyanpink,sections,evestaz_Qc,x,y,z);
    title('Input spike Qc',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    subplot(1,2,2)
    Murat_image3D_2panels(X,Y,Z,spike_outputQc,...
        cyanpink,sections,evestaz_Qc,x,y,z);
    title('Output spike Qc',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    saveas(Qc_spike,fullfile(FPath, FLabel, storeFolder, FNameQcSpike));
    close(Qc_spike)
    %%
    % Spike Q: Input and Output
    FNameQSpike                     =   ['Q-Spike_' fcName '_Hz'];
    Q_spike                         =   figure('Name',FNameQSpike,...
        'NumberTitle','off','Position',[20,400,2000,1000]);
    subplot(1,2,1)
    Murat_image3D_2panels(X,Y,Z,spike_inputQ,...
        purpleorange,sections,evestaz_Q,x,y,z);
    title('Input spike Q',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    subplot(1,2,2)
    Murat_image3D_2panels(X,Y,Z,spike_outputQ,...
        purpleorange,sections,evestaz_Q,x,y,z);
    title('Output spike Q',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    saveas(Q_spike,fullfile(FPath, FLabel, storeFolder, FNameQSpike));
    close(Q_spike)
    
    %% PARAMETER PLOT
    % The final figure is the parameter plot separation.
    storeFolder                     =   'Results';
    %%
    % First Qc and Peak delay are separated in 4 quadrants.
    [param_plot,par,para_map]       =...
        Murat_imageParameters(x,y,z,modv_pd_k,modv_Qc_k);
    xlabel('Qc','FontSize',sizeTitle,'FontWeight','bold','Color','k')
    ylabel('Log. peak delay','FontSize',sizeTitle,...
        'FontWeight','bold','Color','k')
    title('Parameter space plot',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    axis square
    FName_Parameters                =...
        ['Parameter_space_variations_' fcName '_Hz'];
    saveas(param_plot,...
        fullfile(FPath, FLabel, storeFolder, FName_Parameters), fformat);
    close(param_plot)

    %%
    % The second part produces the spatial plot, setting each node to the
    % corresponding color.
    for n=1:length(par(:,1))
        locate_para                 =   para_map(:,1) == par(n,1) &...
            para_map(:,2) == par(n,2) & para_map(:,3) == par(n,3);
        para_map(locate_para,4)     =   par(n,4);
    end
    [~,~,~,para_3D]                 =   Murat_fold(x,y,z,para_map(:,4)-5);
    un_X = unique(para_3D);
    fu                              =   find(un_X==-1,1);
    %%
    % The four options are: (1) high for both (red); (2) low for both (green);
    % (3) high for peak delays only (cyan); (4) high for inverse Qc only
    % (orange).
    if isempty(fu)
        
        cma_para                    =...
            [0.7 0.7 0.7;  0 0.8 0; 0 0.6 1; 1 0.6 0];
        HTick                       =...
            {'Average','Low Scattering Low Absorption',...
            'High Scattering','High Absorption'};
        
    else
        
        cma_para                    =...
            [0.7 0.7 0.7;  0 0.8 0; 0 0.6 1; 1 0.6 0; 1 0 0];
        HTick                       =...
            {'Average','Low Scattering Low Absorption',...
            'High Scattering','High Absorption',...
            'High Scattering and Absorption'};
        
    end
    FName_PMap                      =...
        ['Parameter-Map_' fcName '_Hz'];
    ParaMap                         =...
        Murat_image3D(X,Y,Z,para_3D,...
        'bone',sections,evestaz_Qc,x,y,z,FName_Parameters);
    colormap(cma_para)
    colorbar('Ticks',un_X,'TickLabels',HTick);
    title('Parameter separation map',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    saveas(ParaMap,fullfile(FPath, FLabel, storeFolder, FName_PMap));
    close(ParaMap)
    
    %% SAVE all results as VTK for visualization in PARAVIEW
    storeFolder                     =   'VTK';
    %%
    % Converting Lon/Lat to km for paraview visualization with ndgrid
    x_origin                        =   x - origin(1);
    y_origin                        =   y - origin(2);
    UTM_WE                          =   WE_origin + deg2km(x_origin)*1000;
    UTM_SN                          =   SN_origin + deg2km(y_origin)*1000;
    [X_UTM,Y_UTM,~]                 =   ndgrid(UTM_WE,UTM_SN,z);
    %%
    % Writes the four models to vtk
    vtkwrite(fullfile(FPath, FLabel,storeFolder,[FName_PDMap '.vtk']),...
        'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Peak_delay',mPD)
    vtkwrite(fullfile(FPath, FLabel,storeFolder,[FName_QcMap '.vtk']),...
        'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Qc',mQc)
    vtkwrite(fullfile(FPath, FLabel,storeFolder,[FName_QMap '.vtk']),...
        'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Q',mQ)
    %%
    % Write the input-output checkerboard
    vtkwrite(fullfile(FPath, FLabel,storeFolder,[FNameQcCheck '.vtk']),...
        'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Check_Qc',...
        check_outputQc)
    vtkwrite(fullfile(FPath, FLabel,storeFolder,[FNameQCheck '.vtk']),...
        'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Check_Q',...
        check_outputQ)
    %%
    % Writes the input-output spikes
    vtkwrite(fullfile(FPath, FLabel,storeFolder,[FNameQSpike '.vtk']),...
        'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Spike_Qc',...
        spike_outputQc)
    vtkwrite(fullfile(FPath, FLabel,storeFolder,[FNameQSpike '.vtk']),...
        'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Spike_Q',...
        spike_outputQ)
end

%%
% Also showing the velocity model in case it is available
if Murat.input.availableVelocity == 1
    FName_Vimage                =   'Velocity_model';
    Vimage                      =   Murat_image3D(X,Y,Z,plotV,...
        inferno,sections,evestaz_Q,x,y,z,FName_Vimage);
    title('Velocity Model',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    saveas(Vimage,fullfile(FPath, FLabel, storeFolder, FName_Vimage));
    close(Vimage)
end
vtkwrite(fullfile(FPath, FLabel,storeFolder,'Velocity_model.vtk'),...
    'structured_grid',X_UTM,Y_UTM,Z1,'scalars','V',plotV)
vtkwrite(fullfile(FPath, FLabel,storeFolder,'Input_checkerboards.vtk'),...
    'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Input_check',...
    check_inputQc)
vtkwrite(fullfile(FPath, FLabel,storeFolder,'Input_spikes.vtk'),...
    'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Input_spikes',...
    spike_inputQc)

