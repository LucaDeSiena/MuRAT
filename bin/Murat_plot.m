%% MURAT_PLOT Creates files for visualization in Matlab and Paraview
function Murat                      =   Murat_plot(Murat)
%%
% Importing all the necessary inputs and data for plotting
FLabel                              =   Murat.input.label;
origin                              =   Murat.input.origin;
ending                              =   Murat.input.end;
x                                   =   Murat.input.x;
y                                   =   Murat.input.y;
z                                   =   Murat.input.z;
fformat                             =   Murat.input.format;
sections                            =   Murat.input.sections;
plotV                               =   Murat.input.modvPlot;
cf                                  =   Murat.input.centralFrequency;
vS                                  =   Murat.input.averageVelocityS;
tWm                                 =   Murat.input.codaWindow;
kT                                  =   Murat.input.tresholdNoise;
B0                                  =   Murat.input.albedo;
Le1                                 =   Murat.input.extinctionLength;
QcM                                 =   Murat.input.QcMeasurement;
sped                                =   Murat.input.spectralDecay;
lapseTimeMethod                     =   Murat.input.lapseTimeMethod;

modvQc                              =   Murat.input.modv;
stepgX                              =   (modvQc(2,1) - modvQc(1,1))/2;
stepgY                              =   (modvQc(2,2) - modvQc(1,2))/2;
stepgZ                              =   (modvQc(2,3) - modvQc(1,3))/2;
modvQc(:,1)                         =   modvQc(:,1) + stepgX;
modvQc(:,2)                         =   modvQc(:,2) + stepgY;
modvQc(:,3)                         =   modvQc(:,3) + stepgZ;

Qm                                  =   Murat.data.inverseQc;
time0                               =   Murat.data.travelTime;
retainPeakDelay                     =   Murat.data.retainPeakDelay;
retainQc                            =   Murat.data.retainQc;
retainQ                             =   Murat.data.retainQ;
ray_crosses_pd                      =   Murat.data.raysPeakDelay;
ray_crosses_Qc                      =   Murat.data.raysQc;
ray_crosses_Q                       =   Murat.data.raysQ;
fitrobust                           =   Murat.data.fitrobust;
peakData                            =   Murat.data.peakDelay;
luntot                              =   Murat.data.totalLengthRay;
rma                                 =   Murat.data.raysPlot;
modv_pd                             =   Murat.data.modvPeakDelay;
modv_Qc                             =   Murat.data.modvQc;
modv_Q                              =   Murat.data.modvQ;
evestazDegrees                      =   Murat.data.locationsDeg;
energyRatio                         =   Murat.data.energyRatioBodyCoda;
codaNoiseRatio                      =   Murat.data.energyRatioCodaNoise;
Ac_i                                =   Murat.data.inversionMatrixQc;
RZZ                                 =   Murat.data.uncertaintyQc;
A_i                                 =   Murat.data.inversionMatrixQ;
residualQc                          =   Murat.data.residualQc;
residualQ                           =   Murat.data.residualQ;
locationM                           =   Murat.data.locationsM;
tCoda                               =   Murat.data.tCoda;
rapsp                               =   Murat.data.energyRatioBodyCoda;

FPath                               =   './';
sizeTitle                           =   18;
lMF                                 =   size(ray_crosses_pd);
sections(3)                         =   sections(3)/1000;

cyanpink                            =   colMapGen([1,0,1],[0,1,1],256);
purpleorange                        =...
    colMapGen([0.5 0 0.5],[0.91 0.41 0.17],256);

evestaz                             =...
    [evestazDegrees(:,1:2) -evestazDegrees(:,3)/1000 ...
    evestazDegrees(:,4:5) evestazDegrees(:,6)/1000];

%% PLOTS - coverage and sensitivity
% Murat_plot starts plotting the ray distribution if asked by the user.
% It stores  the files in the corresponding folder.
averageQcFrequency                  =   zeros(2,lMF(2));
for k = 1:lMF(2)
    storeFolder                     =   'Rays_Checks';
    cf_k                            =   cf(k);
    fcName                          =   num2str(cf_k);
    if find(fcName == '.')
        fcName(fcName == '.')       =   '_';
    end
    rtpdk                           =   retainPeakDelay(:,k);
    rtQk                            =   retainQ(:,k);
    rcQk                            =   ray_crosses_Q(:,k);
    rtQck                           =   retainQc(:,k);
    rcQck                           =   ray_crosses_Qc(:,k);
    %%
    % The rays are visualized for different techniques, starting with the peak delay
    FName_peakDelay                 =   ['Rays_PeakDelay_' fcName '_Hz'];
    rma_pd                          =   rma(:,2:4,rtpdk)/1000;
    evestaz_pd                      =   evestaz(rtpdk,:);
    rays_peakDelay                  =   Murat_imageRays(rma_pd,origin,...
        ending,evestaz_pd,x,y,z,FName_peakDelay);
    saveas(rays_peakDelay,fullfile(FPath,FLabel,storeFolder,...
        FName_peakDelay));
    close(rays_peakDelay)
    %%
    % The next figure shows the rays for the total attenuation (Q)
    FName_Q                         =   ['Rays_Q_' fcName '_Hz'];
    
    rma_Q                           =   rma(:,2:4,rtQk)/1000;
    evestaz_Q                       =   evestaz(rtQk,:);
    rays_Q                          =...
        Murat_imageRays(rma_Q,origin,ending,evestaz_Q,...
        x,y,z,FName_Q);
    saveas(rays_Q,fullfile(FPath, FLabel, storeFolder, FName_Q));
    close(rays_Q)
    %%
    % The next sensitivity to check is the one for coda attenuation. The code creates
    % a figure that shows sections in the sensitivity kernels. The left panel shows
    % the sensitivity kernel in the full space while the second shows the normalized
    % kernel in the inversion grid.
    FName_Qc                        =   ['Kernel_Qc' fcName '_Hz'];
    kernels                         =   figure('Name',FName_Qc,...
        'NumberTitle','off','Position',[20,400,1200,1000]);
    
    % Calculates kernels
    [K_grid, r_grid]                =...
        Murat_kernels(tCoda(1)+tWm/2,locationM(1,1:3),locationM(1,4:6),...
        modvQc,vS,kT,B0,Le1,lapseTimeMethod);
    
    Murat_codaMatrix(modvQc,K_grid,r_grid,1,origin,sections);
    
    saveas(kernels,fullfile(FPath, FLabel, storeFolder, FName_Qc));
    close(kernels)
    
    %% Plot - Checks
    % In this section Murat_plot makes checks on the three parameters.
    % These plots are always visualised. They check that:
    % (1) Qc is constant with ray length - also computes weighted average;
    % (2) peak delays increase with travel time;
    % (3) amplitude ratios decay with hypocentral distance.
    % These plots are used to select measurements and understand how well
    % they follow the assumptions.
    Qm_k                            =   Qm(rtQck,k);
    RZZ_k                           =   RZZ(rtQck,k);
    residualQc_k                    =   residualQc(k);
    luntot_Qc                       =   luntot(rtQck)/1000;
    Ac                              =   Ac_i(rtQck,rcQck);

    averageQcFrequency(1,k)         =   sum(RZZ_k.*Qm_k)/sum(RZZ_k);
    averageQcFrequency(2,k)         =   std(Qm_k);
    
    Qc_title                        =   ['Qc check ' fcName ' Hz'];
    Qc_analysis                     =   Murat_imageCheckQc(Qm_k,RZZ_k,...
        residualQc_k,luntot_Qc,Ac,sizeTitle,Qc_title,QcM);
    
    FName                           =   ['Qc_analysis_' fcName '_Hz'];
    saveas(Qc_analysis, fullfile(FPath,FLabel,storeFolder,FName),fformat);
    %%
    % Then it shows the peak delay relative to the travel time.
    peakData_k                      =   peakData(rtpdk,k);
    fitrobust_k                     =   fitrobust(:,k);
    time0PD                         =   time0(rtpdk);
    
    pd_title                        =   ['Peak Delay check ' fcName ' Hz'];
    pd_analysis                     =   Murat_imageCheckPeakDelay(...
    time0PD,fitrobust_k,peakData_k,sizeTitle,pd_title);

    FName                           =   ['PD_analysis_' fcName '_Hz'];
    saveas(pd_analysis, fullfile(FPath,FLabel,storeFolder,FName), fformat);
    %%
    % Then it plots first the logarithm of the energy ratio versus travel
    % time.
    energyRatio_k                   =   energyRatio(rtQk,k);
    residualQ_k                     =   residualQ(k);
    Edirect_k                       =...
        energyRatio_k./codaNoiseRatio(rtQk,k);
    A_k                             =   A_i(rtQk,rcQk);
    luntot_k                        =   luntot(rtQk);
    time0_k                         =   time0(rtQk);
    rapsp_k                         =   rapsp(rtQk,k);    
    tCm                             =   tCoda(rtQk,k);
    Q_k                             =   Qm(rtQk,k);
    
    CN_title                        =...
        ['Coda Normalization check ' fcName ' Hz'];   
    
    [d1, ~,spreadAverageQ, equationQ]...
                                    =...
        Murat_lsqlinQmean(tCm,tWm,Q_k,cf_k,sped,luntot_k,time0_k,rapsp_k);
    
    CN_analysis                     =...
        Murat_imageCheckCN(equationQ,residualQ_k,d1,spreadAverageQ,...
        luntot_k,time0_k,energyRatio_k,A_k,Edirect_k,CN_title);
    
    FName                           =   ['CN_analysis_' fcName '_Hz'];
    saveas(CN_analysis, fullfile(FPath,FLabel,storeFolder,FName), fformat);
    
    %% PLOT - RESULTS
    % Set up matrices. The points are set to the upper SW vertices to
    % work with the function "slice". All stored in the sub-folder.
    storeFolder                     =   'Results';
    modv_pd_k                       =   modv_pd(:,:,k);
    modv_Qc_k                       =   modv_Qc(:,:,k);
    modv_Q_k                        =   modv_Q(:,:,k);
    [X,Y,Z1,mPD]                    =   Murat_fold(x,y,z,modv_pd_k(:,4));
    [~,~,~,mQc]                     =   Murat_fold(x,y,z,modv_Qc_k(:,4));
    [~,~,~,mQ]                      =   Murat_fold(x,y,z,modv_Q_k(:,4));
    Z                               =   Z1/1000;
    evestaz_Qc                      =   evestaz(rtQck,:);
    %%
    % Peak delays results
    FName_PDMap                     =   ['Peak-Delay-3D_' fcName '_Hz'];
    peakDelaymap                    =   Murat_image3D(X,Y,Z,mPD,...
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
    %%
    % First Qc and Peak delay are separated in 4 quadrants.
    % The second part produces the spatial plot, setting each node to the
    %   corresponding color. The four options are: (1) high for both (red);
    %   (2) low for both (green); (3) high for peak delays only (cyan);
    %   (4) high for inverse Qc only (orange).    
    storeFolder                     =   'Results';
    
    %%
    % Define all the parameters for imaging    
    FName_Parameters                =...
        ['Parameter_space_variations_' fcName '_Hz'];
    [param_plot,par,para_map]       =...
        Murat_imageParameters(x,y,z,modv_pd_k,modv_Qc_k,sizeTitle);
    
    saveas(param_plot,...
        fullfile(FPath,FLabel,storeFolder,FName_Parameters), fformat);
    close(param_plot)
    %%
    % Imaging the parameters in 3D
    FName_PMap                      =   ['Parameter-Map_' fcName '_Hz'];
    ParaMap                         = ...
        Murat_imageParametersMaps(par,para_map,x,y,z,X,Y,Z,evestaz_Qc,...
        sections,sizeTitle,FName_Parameters);
    saveas(ParaMap,fullfile(FPath, FLabel, storeFolder, FName_PMap));
    close(ParaMap)
    %% SAVE all results as VTK for visualization in PARAVIEW
    % Converting Lon/Lat to km for paraview visualization with ndgrid
    storeFolder                     =   'VTK';
    [WE_origin, SN_origin]          =   deg2utm(origin(2),origin(1));
    x_origin                        =   x - origin(2);
    y_origin                        =   y - origin(1);
    UTM_WE                          =   WE_origin + deg2km(x_origin)*1000;
    UTM_SN                          =   SN_origin + deg2km(y_origin)*1000;
    [X_UTM,Y_UTM,~]                 =   meshgrid(UTM_WE,UTM_SN,z);
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
% Also showing the velocity model in case it is available, only once
vtkwrite(fullfile(FPath, FLabel,storeFolder,'Velocity_model.vtk'),...
    'structured_grid',X_UTM,Y_UTM,Z1,'scalars','V',plotV)
vtkwrite(fullfile(FPath, FLabel,storeFolder,'Input_checkerboards.vtk'),...
    'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Input_check',...
    check_inputQc)
vtkwrite(fullfile(FPath, FLabel,storeFolder,'Input_spikes.vtk'),...
    'structured_grid',X_UTM,Y_UTM,Z1,'scalars','Input_spikes',...
    spike_inputQc)
%%
% Final figure is the velocity model and relation between Qc and frequency
storeFolder                     =   'Results';

if Murat.input.availableVelocity == 1
    FName_Vimage                    =   'Velocity_model';
    Vimage                          =   Murat_image3D(X,Y,Z,plotV,...
        inferno,sections,evestaz_Q,x,y,z,FName_Vimage);
    title('Velocity Model',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    saveas(Vimage,fullfile(FPath, FLabel, storeFolder, FName_Vimage));
    close(Vimage)
end

Murat.data.averageQcFrequency       =   averageQcFrequency;
Qcf_title                           =   'Qc vs Frequency';
QcFrequency                         =   Murat_imageQcFrequency(cf,...
    averageQcFrequency,sizeTitle,Qcf_title);

FName                               =   'Qc_vs_frequency';
saveas(QcFrequency, fullfile(FPath,FLabel,storeFolder,FName), fformat);
close(QcFrequency)
