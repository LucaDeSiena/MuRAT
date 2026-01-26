%% MURAT_PLOT Creates files for visualization in Matlab
function Murat      =   Murat_plot(Murat)

% Refactored input unpacking and preallocation
s               =   Murat.input;

% Importing all the necessary inputs and data for plotting
FLabel          =   s.label;
origin          =   s.origin;
ending          =   s.end;
x               =   s.x;
y               =   s.y;
z               =   s.z;
sections        =   s.sections;
cf              =   s.centralFrequency;
tWm             =   s.codaWindow;
D               =   s.DiffConstant;
PlotRays        =   s.PlotRays;
PlotTests       =   s.PlotTests;
PlotResults     =   s.PlotResults;
PlotCheckers    =   s.PlotCheckers;
PlotSpikes      =   s.PlotSpikes;
PlotParameters  =   s.PlotParameters;
CoordGrid       =   s.DDcoordinates;

luntot          =   Murat.rays.totalLengthRay;
rma             =   Murat.rays.raysPlot;
time0           =   Murat.rays.travelTime;
evstDegrees     =   Murat.rays.locationsDeg;
retainPeakDelay =   Murat.PD.retainPeakDelay;
rayChPD         =   Murat.PD.raysPeakDelay;
fitrobust       =   Murat.PD.fitrobust;
peakData        =   Murat.PD.peakDelay;
modv_pd         =   Murat.PD.modvPeakDelay;
Qm              =   Murat.Qc.inverseQc;
retainQc        =   Murat.Qc.retainQc;
Ac_i            =   Murat.Qc.inversionMatrixQc;
RZZ             =   Murat.Qc.uncertaintyQc;
residualQc      =   Murat.Qc.residualRedQc;
tCoda           =   Murat.Qc.tCoda;
ray_crosses_Qc  =   Murat.Qc.raysQc;
modv_Qc         =   Murat.Qc.modvQc;
retainQ         =   Murat.Q.retainQ;
modv_Q          =   Murat.Q.modvQ;
energyRatio     =   Murat.Q.energyRatioBodyCoda;
codaNoiseRatio  =   Murat.Q.energyRatioCodaNoise;

sTitle          =   18;

lMF             =   size(rayChPD);
sections(3)     =   sections(3)/1000;

%% PLOTS - coverage and sensitivity
evst            =  [evstDegrees(:,1:2) -evstDegrees(:,3)/1000 ...
    evstDegrees(:,4:5) evstDegrees(:,6)/1000];

avQcFreq        =   zeros(2,lMF(2));

% small helper to build and save path once
makePath = @(varargin) fullfile('./', FLabel, varargin{:});

for k = 1:lMF(2)
    %%
    D_k         =   D(k);
    cf_k        =   cf(k);
    fcName      = strrep(num2str(cf_k),'.','_');
    rtpdk       =   retainPeakDelay(:,k);
    rtQk        =   retainQ(:,k);
    rtQck       =   retainQc(:,k);
    rcQck       =   ray_crosses_Qc(:,k);
    
    if PlotRays == 1
    % Murat_plot starts plotting the ray distribution if asked by the user.
    % It stores  the files in the corresponding folder.
    storeFolder =   'Rays';
    
    % Peak Delay rays
    FName_peakDelay = ['Rays_PeakDelay_' fcName '_Hz'];
    rma_pd      = rma(:,2:4,rtpdk)/1000;
    evst_pd     =   evst(rtpdk,:);
    rays_peakDelay  =   Murat_imageRays(rma_pd,origin,ending,evst_pd,...
        x,y,z,FName_peakDelay);
    saveFigureAsImage(makePath(storeFolder, FName_peakDelay));   
    close(rays_peakDelay)
    
    %%
    % The next figure shows the rays for the total attenuation (Q)
    FName_Q     =   ['Rays_Q_' fcName '_Hz'];
    rma_Q       =   rma(:,2:4,rtQk)/1000;
    evst_Q      =   evst(rtQk,:);
    rays_Q      =   Murat_imageRays(rma_Q,origin,ending,evst_Q,x,y,z,...
        FName_Q);
    saveFigureAsImage(makePath(storeFolder, FName_Q));
    close(rays_Q)

    end
    % Tests
    % In this section Murat_plot makes checks on the three parameters.
    % These plots are always visualised. They check that:
    % (1) Qc is constant with ray length - also computes weighted average;
    % (2) peak delays increase with travel time;
    % (3) amplitude ratios decay with hypocentral distance.
    % These plots are used to select measurements and understand how well
    % they follow the assumptions.

    if PlotTests == 1
    
    % Qc test
    storeFolder =   fullfile('Tests','Qc');
    Qm_k        =   Qm(rtQck,k);
    RZZ_k       =   RZZ(rtQck,k);
    residualQc_k=   residualQc(k);
    luntot_Qc   =   luntot(rtQck)/1000;
    Ac          =   Ac_i(rtQck,rcQck);

    avQcFreq(1,k)   =   sum(RZZ_k.*Qm_k)/sum(RZZ_k);
    avQcFreq(2,k)   =   std(Qm_k);

    Qc_title    =   ['Qc check ' fcName ' Hz'];
    Qc_analysis =   Murat_imageCheckQc(Qm_k,RZZ_k,residualQc_k,...
        luntot_Qc,Ac,sTitle,Qc_title);
    p           =   makePath(storeFolder, ['Qc_analysis_' fcName '_Hz']);
    saveFigureAsImage(p);
    savefig(Qc_analysis, [p '.fig']);
    close(Qc_analysis);
    
    % Peak delay test
    storeFolder =   fullfile('Tests','PeakDelay');
    peakData_k  =   peakData(rtpdk,k);
    fitrobust_k =   fitrobust(:,k);
    time0PD     =   time0(rtpdk);

    pd_title    =   ['Peak Delay check ' fcName ' Hz'];
    pd_analysis =   Murat_imageCheckPeakDelay(time0PD,fitrobust_k,...
        peakData_k,sTitle,pd_title);
    p           = makePath(storeFolder, ['PD_analysis_' fcName '_Hz']);
    saveFigureAsImage(p);
    savefig(pd_analysis, [p '.fig']);
    close(pd_analysis);

    % Coda normalization test
    storeFolder =   fullfile('Tests','Q');
    mask        =   rtQk & rtQck;
        
    energyRatio_k   =   energyRatio(mask,k);
    Ed_k        =   energyRatio_k./codaNoiseRatio(mask,k);
    Qc_k        =   Qm(mask,k);
    luntot_k    =   luntot(mask);
    time0_k     =   time0(mask);
    rapsp_k     =   energyRatio(mask,k);
    tCm         =   tCoda(mask,k);
    l           =   luntot_k/1000;

    CNTitle     =   ['Coda Normalization check ' fcName ' Hz'];
    [d0, Q0, ~] =   Murat_lsqlinQmean(tCm,tWm,Qc_k,cf_k,D_k,l,time0_k,...
        rapsp_k);
    CN_analysis =   Murat_imageCheckCN(time0_k,Q0,d0,Ed_k,CNTitle);
    p           =  makePath(storeFolder,['CN_analysis_' fcName '_Hz']);
    saveFigureAsImage(p);
    savefig(CN_analysis, [p '.fig']);
    close(CN_analysis);

    end

    % PLOT - RESULTS
    if PlotResults == 1
    % Set up matrices. The points are set to the upper SW vertices to
    % work with the function "slice". All stored in the sub-folder.
    storeFolder = fullfile('Results','PeakDelay');

    modv_pd_k   =   modv_pd(:,:,k);
    modv_Qc_k   =   modv_Qc(:,:,k);
    modv_Q_k    =   modv_Q(:,:,k);
    [X,Y,Z1,mPD]=   Murat_fold(x,y,z,modv_pd_k(:,4));
    [~,~,~,mQc] =   Murat_fold(x,y,z,modv_Qc_k(:,4));
    [~,~,~,mQ]  =   Murat_fold(x,y,z,modv_Q_k(:,4));
    Z           =   Z1/1000;
    evst_Qc     =   evst(rtQck,:);
    
    % Peak delays results, using interpolation defined by 'divi'.
    divi            =   5;
    FName_PDMap     =   ['Peak-Delay-3D_' fcName '_Hz'];
    peakDelaymap    =   Murat_image3D(X,Y,Z,mPD,redblue,sections,...
        evst_pd,x,y,z,divi,FName_PDMap);
    title('Log. peak-delay variations','FontSize',sTitle,...
        'FontWeight','bold','Color','k');
    Murat_saveFigures(peakDelaymap, makePath(storeFolder, FName_PDMap));

    % Plots peak delays only keeping cells with more than 'factor'% of data
    factor          =   5;
    max_pos         =   max(mPD(mPD >= 0)); 
    max_neg         =   min(mPD(mPD < 0)); 
    thresholdPlus   =   (max_pos / 100) * factor; 
    thresholdMinus  =   (max_neg / 100) * factor;  
    keepBinsPositive= mPD > thresholdPlus;  
    keepBinsNegative= mPD < thresholdMinus; 
    keep_bins       = keepBinsPositive | keepBinsNegative; 
    
    mPDRed          =   mPD.*keep_bins;

    FName_PDMap     =...
        ['Peak-Delay-3D_' fcName '_Hz_',num2str(factor),'_perc'];
    [peakDelaymapRed,pd_inter,~,~,~,Xi,Yi,Zi] =   Murat_image3D(X,Y,Z,...
        mPDRed,redblue,sections,evst_pd,x,y,z,divi,FName_PDMap);
    title('Peak-delay variations','FontSize',sTitle,...
        'FontWeight','bold','Color','k');
    Murat_saveFigures(peakDelaymapRed, makePath(storeFolder, FName_PDMap));

    % interpolated for the parameter map
    interp_modv_pd_k=   Murat_unfold(Xi,Yi,Zi,pd_inter);

    % Qc results
    storeFolder     =   fullfile('Results','Qc');
    FName_QcMap     =   ['Qc-3D_' fcName '_Hz'];
    [Qcmap,qc_inter,xi,yi,zi,Xi,Yi,Zi] =   Murat_image3D(X,Y,Z,mQc,...
        turbo,sections,evst_Qc,x,y,z,divi,FName_QcMap);
    title('Coda attenuation','FontSize',sTitle,'FontWeight','bold',...
        'Color','k');
    
    Murat_saveFigures(Qcmap, makePath(storeFolder, FName_QcMap));
        
    % interpolated for the parameter map
    interp_modvQc_k =   Murat_unfold(Xi,Yi,Zi,qc_inter);

    %%
    % Q results
    storeFolder     =   fullfile('Results','Q');
    FName_QMap      =   ['Q-3D_' fcName '_Hz'];
    Qmap            =   Murat_image3D(X,Y,Z,mQ,hot,sections,evst_Q,...
        x,y,z,divi,FName_QMap);
    title('Total attenuation variations','FontSize',sTitle,...
        'FontWeight','bold','Color','k');
    Murat_saveFigures(Qmap, makePath(storeFolder, FName_QMap));
    
    end

    %CHECKERBOARDS and HIT MAPS
    % In this section Murat_plot shows the hit map for peak delays and
    % checkerboard tests for Q and Qc.
    if PlotCheckers == 1
    
    storeFolder =   fullfile('Checkerboard','PD');
    FName_PDC   =   ['PD-HitMap_' fcName '_Hz'];
    rPD         =   modv_pd_k(:,5);
    Murat_hitmap(CoordGrid,rPD,makePath(storeFolder, FName_PDC));

    storeFolder =   fullfile('Checkerboard','Qc');
        
    [~,~,~,check_inputQc]   =   Murat_fold(x,y,z,modv_Qc_k(:,6));
    [~,~,~,check_outputQc]  =   Murat_fold(x,y,z,modv_Qc_k(:,7));
    [~,~,~,check_inputQ]    =   Murat_fold(x,y,z,modv_Qc_k(:,6));
    [~,~,~,check_outputQ]   =   Murat_fold(x,y,z,modv_Q_k(:,7));
    
    FName_QcC   =   ['Qc-Checkerboard_' fcName '_Hz'];
    Qc_check    =   myfig(FName_QcC);

    subplot(1,2,1)
    Murat_image3D_2panels(X,Y,Z,check_inputQc,'bone',sections,...
        evst_Qc,x,y,z);
    title('Input checkerboard Qc','FontSize',sTitle,...
        'FontWeight','bold','Color','k');

    subplot(1,2,2)
    Murat_image3D_2panels(X,Y,Z,check_outputQc,...
        'bone',sections,evst_Qc,x,y,z);
    title('Output checkerboard Qc','FontSize',sTitle,...
        'FontWeight','bold','Color','k');

    Murat_saveFigures_2panels(Qc_check,makePath(storeFolder,FName_QcC));

    %Checkerboard Q: Input and Output
    storeFolder     =   fullfile('Checkerboard','Q');
    FName_QCheck    =   ['Q-Checkerboard_' fcName '_Hz'];
    Q_check         =   myfig(FName_QCheck);

    subplot(1,2,1)
    Murat_image3D_2panels(X,Y,Z,check_inputQ,...
        'bone',sections,evst_Q,x,y,z);
    title('Input checkerboard Q','FontSize',sTitle,...
        'FontWeight','bold','Color','k');

    subplot(1,2,2)
    Murat_image3D_2panels(X,Y,Z,check_outputQ,...
        'bone',sections,evst_Q,x,y,z);
    title('Output checkerboard Q',...
        'FontSize',sTitle,'FontWeight','bold','Color','k');
    Murat_saveFigures_2panels(Q_check,makePath(storeFolder, FName_QCheck));
    
    end
    %SPIKES
    
    % In this section Murat_plot shows input and output of the spike tests
    % for Q and Qc.
    
    if PlotSpikes == 1

    [~,~,~,spike_inQc]   =   Murat_fold(x,y,z,modv_Qc_k(:,8));
    [~,~,~,spike_outQc]  =   Murat_fold(x,y,z,modv_Qc_k(:,9));
    [~,~,~,spike_inQ]    =   Murat_fold(x,y,z,modv_Qc_k(:,8));
    [~,~,~,spike_outQ]   =   Murat_fold(x,y,z,modv_Q_k(:,9));

    % Spike Qc: Input and Output
    storeFolder     =   fullfile('Spike','Qc');
    FNameQcSpike    =   ['Qc-Spike_' fcName '_Hz'];
    Qc_spike        =   myfig(FNameQcSpike);

    subplot(1,2,1)
    Murat_image3D_2panels(X,Y,Z,spike_inQc,winter,sections,evst_Qc,x,y,z);
    title('Input spike Qc','FontSize',sTitle,...
        'FontWeight','bold','Color','k');
    subplot(1,2,2)
    Murat_image3D_2panels(X,Y,Z,spike_outQc,...
        winter,sections,evst_Qc,x,y,z);
    title('Output spike Qc','FontSize',sTitle,...
        'FontWeight','bold','Color','k');

    Murat_saveFigures_2panels(Qc_spike,makePath(storeFolder,FNameQcSpike));

    % Spike Q: Input and Output
    storeFolder     =   fullfile('Spike','Q');
    FNameQSpike     =   ['Q-Spike_' fcName '_Hz'];
    Q_spike         =   myfig(FNameQSpike);

    subplot(1,2,1)
    Murat_image3D_2panels(X,Y,Z,spike_inQ,hot,sections,evst_Q,x,y,z);
    title('Input spike Q','FontSize',sTitle,...
        'FontWeight','bold','Color','k');
    subplot(1,2,2)
    Murat_image3D_2panels(X,Y,Z,spike_outQ,...
        hot,sections,evst_Q,x,y,z);
    title('Output spike Q','FontSize',sTitle,...
        'FontWeight','bold','Color','k');

    Murat_saveFigures_2panels(Q_spike,makePath(storeFolder,FNameQSpike));
    end
    
    %% PARAMETER PLOT
    if PlotParameters == 0 %%KILLING THIS FOR NOW
    % The final figure is the parameter plot separation.
    % First Qc and Peak delay are separated in 4 quadrants.
    % The second part produces the spatial plot, setting each node to the
    %   corresponding color. The four options are: (1) high for both (red);
    %   (2) low for both (green); (3) high for peak delays only (cyan);
    %   (4) high for inverse Qc only (orange).
    storeFolder     =   fullfile('Results','Parameter');
    
    % Define all the parameters for imaging
    FName_Param     =   ['Parameter_space_variations_' fcName '_Hz'];
    [param_plot,~,~]=   Murat_imageParameters(x,y,z,modv_pd_k,...
        modv_Qc_k,sTitle);
    savefig(param_plot,makePath(storeFolder,FName_Param));
    close(param_plot)

    % Use interpolated peakdelay and Qc
    zi              =   (zi*1000)';
    [~,par_inter,para_map_inter]    =   Murat_imageParameters(xi',yi',...
        zi,interp_modv_pd_k,interp_modvQc_k,sTitle);

    %%
    % Imaging the parameters in 3D
    FName_PMap      =   ['Parameter-Map_' fcName '_Hz'];
    [ParaMap,para_map]              =...
        Murat_imageParametersMaps(par_inter,para_map_inter,xi',yi',zi,...
        Xi,Yi,Zi,evst_Qc,sections,sTitle,FName_PMap);
    Murat_saveFigures(ParaMap,makePath(storeFolder, FName_PMap));
    
    FName           =   ['parameterMap_' fcName '_Degrees_Hz.txt'];
    writematrix(para_map,makePath('TXT',FName));

    end
end

% Final figure is the Qc vs frequency relation
Murat.Qc.averageQcFrequency = avQcFreq;
Qcf_title = 'Qc vs Frequency';
Murat_imageQcFrequency(cf, avQcFreq, sTitle, Qcf_title);
saveFigureAsImage(makePath('Results', 'Qc_vs_frequency'));
close all;

end