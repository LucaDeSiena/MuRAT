%% MURAT_PLOT Creates files for visualization in Matlab
function Murat                      =   Murat_plot(Murat)
%%
% Importing all the necessary inputs and data for plotting
FLabel                              =   Murat.input.label;
origin                              =   Murat.input.origin;
ending                              =   Murat.input.end;
x                                   =   Murat.input.x;
y                                   =   Murat.input.y;
z                                   =   Murat.input.z;
sections                            =   Murat.input.sections;
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

%% PLOTS - coverage and sensitivity
% Declustering is done before any frequency analysis, here we show the 2D
% rays before and after
if Murat.input.declustering > 0
    locDegOriginal                  =   Murat.data.locDegOriginal;
    FName_Cluster                   =   'Clustering';
    clustering                      =   Murat_imageDeclustering(...
        locDegOriginal,evestazDegrees,origin,ending,FName_Cluster);
    storeFolder                     =   'Tests';
    pathFolder                      =...
        fullfile(FPath,FLabel,storeFolder,FName_Cluster);
    saveas(clustering,pathFolder,'tif');
    close(clustering)
end

evestaz                             =...
    [evestazDegrees(:,1:2) -evestazDegrees(:,3)/1000 ...
    evestazDegrees(:,4:5) evestazDegrees(:,6)/1000];

averageQcFrequency                  =   zeros(2,lMF(2));

for k = 1:lMF(2)
    % Murat_plot starts plotting the ray distribution if asked by the user.
    % It stores  the files in the corresponding folder.
    storeFolder                     =   'Rays';
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
    pathFolder                      =...
        fullfile(FPath,FLabel,storeFolder,FName_peakDelay);
    saveas(rays_peakDelay,pathFolder,'tif');
    close(rays_peakDelay)
    
    %%
    % The next figure shows the rays for the total attenuation (Q)
    FName_Q                         =   ['Rays_Q_' fcName '_Hz'];
    rma_Q                           =   rma(:,2:4,rtQk)/1000;
    evestaz_Q                       =   evestaz(rtQk,:);
    rays_Q                          =...
        Murat_imageRays(rma_Q,origin,ending,evestaz_Q,x,y,z,FName_Q);
    pathFolder                      =...
        fullfile(FPath, FLabel, storeFolder, FName_Q);
    saveas(rays_Q,pathFolder,'tif');
    close(rays_Q)
    
    %%
    % The next figure checks the sensitivity of coda attenuation
    % measurements. The code creates figures that show sections in the
    % sensitivity kernels. The left panel shows the sensitivity kernel in
    % the full space while the rigth panel shows the normalized
    % kernel in the inversion grid.
    storeFolder                     =   'Kernels';
    FName_Qc                        =   ['Kernel_Qc' fcName '_Hz'];
    kernels                         =   figure('Name',FName_Qc,...
        'NumberTitle','off','Position',[20,400,1200,1000],'visible','off');

    % Calculates kernels
    [K_grid, r_grid]                =...
        Murat_kernels(tCoda(1)+tWm/2,locationM(1,1:3),locationM(1,4:6),...
        modvQc,vS,kT,B0,Le1,lapseTimeMethod);

    Murat_codaMatrix(modvQc,K_grid,r_grid,1,origin,sections);
    pathFolder                      =...
        fullfile(FPath, FLabel, storeFolder, FName_Qc);
    Murat_saveFigures_2panels(kernels,pathFolder);
    
    %% Plot - Tests
    % In this section Murat_plot makes checks on the three parameters.
    % These plots are always visualised. They check that:
    % (1) Qc is constant with ray length - also computes weighted average;
    % (2) peak delays increase with travel time;
    % (3) amplitude ratios decay with hypocentral distance.
    % These plots are used to select measurements and understand how well
    % they follow the assumptions.
    storeFolder                     =   'Tests/Qc';
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
    saveas(Qc_analysis, fullfile(FPath,FLabel,storeFolder,...
        ['Qc_analysis_' fcName '_Hz']),'tif');
    saveas(Qc_analysis, fullfile(FPath,FLabel,storeFolder,...
        ['Qc_analysis_' fcName '_Hz']));
    close(Qc_analysis)
    
    %%
    % Then it shows the peak delay relative to the travel time.
    storeFolder                     =   'Tests/PeakDelay';
    peakData_k                      =   peakData(rtpdk,k);
    fitrobust_k                     =   fitrobust(:,k);
    time0PD                         =   time0(rtpdk);

    pd_title                        =   ['Peak Delay check ' fcName ' Hz'];
    pd_analysis                     =   Murat_imageCheckPeakDelay(...
    time0PD,fitrobust_k,peakData_k,sizeTitle,pd_title);
    saveas(pd_analysis, fullfile(FPath,FLabel,storeFolder,...
        ['PD_analysis_' fcName '_Hz']),'tif');
    saveas(pd_analysis, fullfile(FPath,FLabel,storeFolder,...
        ['PD_analysis_' fcName '_Hz']));
    close(pd_analysis)
    
    %%
    % Then it plots first the logarithm of the energy ratio versus travel
    % time.
    storeFolder                     =   'Tests/Q';
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
    saveas(CN_analysis, fullfile(FPath,FLabel,storeFolder,...
        ['CN_analysis_' fcName '_Hz']),'tif');
    saveas(CN_analysis, fullfile(FPath,FLabel,storeFolder,...
        ['CN_analysis_' fcName '_Hz']));
    close(CN_analysis)
    
    %% PLOT - RESULTS
    % Set up matrices. The points are set to the upper SW vertices to
    % work with the function "slice". All stored in the sub-folder.
    modv_pd_k                       =   modv_pd(:,:,k);
    modv_Qc_k                       =   modv_Qc(:,:,k);
    modv_Q_k                        =   modv_Q(:,:,k);
    [X,Y,Z1,mPD]                    =   Murat_fold(x,y,z,modv_pd_k(:,4));
    [~,~,~,PD_cts]                  =   Murat_fold(x,y,z,modv_pd_k(:,5));
    [~,~,~,mQc]                     =   Murat_fold(x,y,z,modv_Qc_k(:,4));
    [~,~,~,mQ]                      =   Murat_fold(x,y,z,modv_Q_k(:,4));
    Z                               =   Z1/1000;
    evestaz_Qc                      =   evestaz(rtQck,:);
    
    %%
    % Peak delays results, using interpolation defined by 'divi'.
    divi                            =   5;
    storeFolder                     =   'Results/PeakDelay';
    FName_PDMap                     =   ['Peak-Delay-3D_' fcName '_Hz'];
    peakDelaymap                    =   Murat_image3D(X,Y,Z,mPD,...
        redblue,sections,evestaz_pd,x,y,z,divi,FName_PDMap);
    title('Log. peak-delay variations',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    pathFolder                      =...
        fullfile(FPath, FLabel, storeFolder, FName_PDMap);
    Murat_saveFigures(peakDelaymap,pathFolder);

    %%
    % Plots peak delays only keeping cells with more than 'factor'% of data
    factor                          =   5;
    keep_bins                       =...
        PD_cts > ((max(PD_cts(:))/100)*factor);
    mPD_red                         =   mPD.*keep_bins;

    FName_PDMap                     =...
        ['Peak-Delay-3D_' fcName '_Hz_',num2str(factor),'_perc'];
    [peakDelaymap_red,pd_inter,~,~,~,Xi,Yi,Zi]...
                                    =   Murat_image3D(X,Y,Z,mPD_red,...
        redblue,sections,evestaz_pd,x,y,z,divi,FName_PDMap);
    title('Peak-delay variations','FontSize',sizeTitle,...
        'FontWeight','bold','Color','k');
    pathFolder                      =...
        fullfile(FPath, FLabel, storeFolder, FName_PDMap);
    Murat_saveFigures(peakDelaymap_red,pathFolder);
    
    % interpolated for the parameter map
    interp_modv_pd_k                =   Murat_unfold(Xi,Yi,Zi,pd_inter);

    %%
    % Qc results
    storeFolder                     =   'Results/Qc';
    FName_QcMap                     =   ['Qc-3D_' fcName '_Hz'];
    [Qcmap,qc_inter,xi,yi,zi,Xi,Yi,Zi]...
                                    =   Murat_image3D(X,Y,Z,mQc,...
        turbo,sections,evestaz_Qc,x,y,z,divi,FName_QcMap);
    title('Coda attenuation',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');
    pathFolder                      =...
        fullfile(FPath,FLabel,storeFolder,FName_QcMap);
    Murat_saveFigures(Qcmap,pathFolder);

    % interpolated for the parameter map
    interp_modv_qc_k                =   Murat_unfold(Xi,Yi,Zi,qc_inter);

    %%
    % Q results
    storeFolder                     =   'Results/Q';
    FName_QMap                      =   ['Q-3D_' fcName '_Hz'];
    Qmap                            =   Murat_image3D(X,Y,Z,mQ,...
        hot,sections,evestaz_Q,x,y,z,divi,FName_QMap);
    title('Total attenuation variations','FontSize',sizeTitle,...
        'FontWeight','bold','Color','k');
    pathFolder                      =...
        fullfile(FPath,FLabel,storeFolder,FName_QMap);
    Murat_saveFigures(Qmap,pathFolder);

    %% PLOT - CHECKERBOARDS
    % In this section Murat_plot shows the checkerboard tests
    % for Q and Qc.
    [~,~,~,check_inputQc]           =   Murat_fold(x,y,z,modv_Qc_k(:,6));
    [~,~,~,check_outputQc]          =   Murat_fold(x,y,z,modv_Qc_k(:,7));
    [~,~,~,check_inputQ]            =   Murat_fold(x,y,z,modv_Qc_k(:,6));
    [~,~,~,check_outputQ]           =   Murat_fold(x,y,z,modv_Q_k(:,7));
    
    %%
    % Checkerboard Qc: Input and Output
    storeFolder                     =   'Checkerboard/Qc';
    FName_QcCheck                   =   ['Qc-Checkerboard_' fcName '_Hz'];
    Qc_check                        =   figure('Name',FName_QcCheck,...
        'NumberTitle','off','Position',[20,400,2000,1000],'visible','off');

    subplot(1,2,1)
    Murat_image3D_2panels(X,Y,Z,check_inputQc,...
        'bone',sections,evestaz_Qc,x,y,z);
    title('Input checkerboard Qc','FontSize',sizeTitle,...
        'FontWeight','bold','Color','k');

    subplot(1,2,2)
    Murat_image3D_2panels(X,Y,Z,check_outputQc,...
        'bone',sections,evestaz_Qc,x,y,z);
    title('Output checkerboard Qc','FontSize',sizeTitle,...
        'FontWeight','bold','Color','k');

    pathFolder                      =...
        fullfile(FPath,FLabel,storeFolder,FName_QcCheck);
    Murat_saveFigures_2panels(Qc_check,pathFolder);

    %%
    %Checkerboard Q: Input and Output
    storeFolder                     =   'Checkerboard/Q';
    FName_QCheck                    =   ['Q-Checkerboard_' fcName '_Hz'];
    Q_check                         =   figure('Name',FName_QCheck,...
        'NumberTitle','off','Position',[20,400,2000,1000],'visible','off');

    subplot(1,2,1)
    Murat_image3D_2panels(X,Y,Z,check_inputQ,...
        'bone',sections,evestaz_Q,x,y,z);
    title('Input checkerboard Q','FontSize',sizeTitle,...
        'FontWeight','bold','Color','k');

    subplot(1,2,2)
    Murat_image3D_2panels(X,Y,Z,check_outputQ,...
        'bone',sections,evestaz_Q,x,y,z);
    title('Output checkerboard Q',...
        'FontSize',sizeTitle,'FontWeight','bold','Color','k');

    pathFolder                      =...
        fullfile(FPath,FLabel,storeFolder,FName_QCheck);
    Murat_saveFigures_2panels(Q_check,pathFolder);

    %% PLOT - SPIKES
    % In this section Murat_plot shows input and output of the spike tests
    % for Q and Qc.
    [~,~,~,spike_inputQc]           =   Murat_fold(x,y,z,modv_Qc_k(:,8));
    [~,~,~,spike_outputQc]          =   Murat_fold(x,y,z,modv_Qc_k(:,9));
    [~,~,~,spike_inputQ]            =   Murat_fold(x,y,z,modv_Qc_k(:,8));
    [~,~,~,spike_outputQ]           =   Murat_fold(x,y,z,modv_Q_k(:,9));
    
    %%
    % Spike Qc: Input and Output
    storeFolder                     =   'Spike/Qc';
    FName_QcSpike                   =   ['Qc-Spike_' fcName '_Hz'];
    Qc_spike                        =   figure('Name',FName_QcSpike,...
        'NumberTitle','off','Position',[20,400,2000,1000],'visible','off');

    subplot(1,2,1)
    Murat_image3D_2panels(X,Y,Z,spike_inputQc,...
        winter,sections,evestaz_Qc,x,y,z);
    title('Input spike Qc','FontSize',sizeTitle,...
        'FontWeight','bold','Color','k');
    subplot(1,2,2)
    Murat_image3D_2panels(X,Y,Z,spike_outputQc,...
        winter,sections,evestaz_Qc,x,y,z);
    title('Output spike Qc','FontSize',sizeTitle,...
        'FontWeight','bold','Color','k');

    pathFolder                      =...
        fullfile(FPath,FLabel,storeFolder,FName_QcSpike);
    Murat_saveFigures_2panels(Qc_spike,pathFolder);

    %%
    % Spike Q: Input and Output
    storeFolder                     =   'Spike/Q';
    FName_QSpike                    =   ['Q-Spike_' fcName '_Hz'];
    Q_spike                         =   figure('Name',FName_QSpike,...
        'NumberTitle','off','Position',[20,400,2000,1000],'visible','off');

    subplot(1,2,1)
    Murat_image3D_2panels(X,Y,Z,spike_inputQ,hot,sections,evestaz_Q,x,y,z);
    title('Input spike Q','FontSize',sizeTitle,...
        'FontWeight','bold','Color','k');
    subplot(1,2,2)
    Murat_image3D_2panels(X,Y,Z,spike_outputQ,...
        hot,sections,evestaz_Q,x,y,z);
    title('Output spike Q','FontSize',sizeTitle,...
        'FontWeight','bold','Color','k');

    pathFolder                      =...
        fullfile(FPath,FLabel,storeFolder,FName_QSpike);
    Murat_saveFigures_2panels(Q_spike,pathFolder);

    %% PARAMETER PLOT
    % The final figure is the parameter plot separation.
    % First Qc and Peak delay are separated in 4 quadrants.
    % The second part produces the spatial plot, setting each node to the
    %   corresponding color. The four options are: (1) high for both (red);
    %   (2) low for both (green); (3) high for peak delays only (cyan);
    %   (4) high for inverse Qc only (orange).
    storeFolder                     =   'Results/Parameter';
    
    %%
    % Define all the parameters for imaging
    FName_Parameters                =...
        ['Parameter_space_variations_' fcName '_Hz'];
    [param_plot,~,~]       =...
        Murat_imageParameters(x,y,z,modv_pd_k,modv_Qc_k,sizeTitle);
    saveas(param_plot,fullfile(FPath,FLabel,storeFolder,FName_Parameters));
    close(param_plot)

    % use interpolated peakdelay and Qc
    zi = (zi*1000)';
    [~,par_inter,para_map_inter]    =   Murat_imageParameters(xi',yi',...
        zi,interp_modv_pd_k,interp_modv_qc_k,sizeTitle);

    %%
    % Imaging the parameters in 3D
    FName_PMap                      =   ['Parameter-Map_' fcName '_Hz'];
    [ParaMap,para_map]              =...
        Murat_imageParametersMaps(par_inter,para_map_inter,xi',yi',zi,...
        Xi,Yi,Zi,evestaz_Qc,sections,sizeTitle,FName_PMap);
    
    pathFolder                      =...
        fullfile(FPath,FLabel,storeFolder,FName_PMap);
    Murat_saveFigures(ParaMap,pathFolder);
    
    FName                           =...
        ['parameterMap_' fcName '_Degrees_Hz.txt'];
    writematrix(para_map,fullfile(FPath, FLabel, 'TXT', FName));
end

%%
% Final figure is the Qc vs frequency relation
storeFolder                         =   'Tests';
Murat.data.averageQcFrequency       =   averageQcFrequency;
Qcf_title                           =   'Qc vs Frequency';
QcFrequency                         =   Murat_imageQcFrequency(cf,...
    averageQcFrequency,sizeTitle,Qcf_title);
FName                               =   'Qc_vs_frequency';
saveas(QcFrequency, fullfile(FPath,FLabel,storeFolder,FName),'tif');
close all
