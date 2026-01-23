%% Peak-delay, Qc and Q TOMOGRAPHIC INVERSIONS
function Murat                  =   Murat_inversion(Murat)
% Importing all the necessary inputs and data for plotting

s                   =   Murat.input;  % local shorthand

% Inputs
FLabel              =   s.label;
tWm                 =   s.codaWindow;
cf                  =   s.centralFrequency;
sizea               =   s.sizeCheck;
modvP               =   s.modvPlot;
spike_o             =   s.spikeLocationOrigin;
spike_e             =   s.spikeLocationEnd;
spike_v             =   s.spikeValue;
x                   =   s.x;
y                   =   s.y;
z                   =   s.z;
nxc                 =   numel(x);
nyc                 =   numel(y);
nzc                 =   numel(z);
nxyzc               =   nxc*nyc*nzc;
inversionMethod     =   s.inversionMethod;
muratHeader         =   s.header;
DDcoordinates       =   s.DDcoordinates;
B0                  =   s.albedo;
Le_1                =   s.iExtinctionLength;
vS                  =   s.averageVelocityS;
iter                =   s.MaximumIterations;
iterStall           =   s.MaximumStallIterations;
dampValueQc         =   s.dampingValueQc;
smoothValueQc       =   s.smoothingValueQc;
dampValueQ          =   s.dampingValueQ;
smoothValueQ        =   s.smoothingValueQ;
plotI               =   s.PlotInversion;

% Matrices / data from Murat struct
Apd_i               =   Murat.PD.inversionMatrixPeakDelay;
Ac_i                =   Murat.Qc.inversionMatrixQc;
A_i                 =   Murat.Q.inversionMatrixQ;
luntot              =   Murat.rays.totalLengthRay;
time0               =   Murat.rays.travelTime;
Qm                  =   Murat.Qc.inverseQc;
RZZ                 =   Murat.Qc.uncertaintyQc;
lpdelta             =   Murat.PD.variationPeakDelay;
rapsp               =   Murat.Q.energyRatioBodyCoda;
retain_pd           =   Murat.PD.retainPeakDelay;
retain_Qc           =   Murat.Qc.retainQc;
retain_Q            =   Murat.Q.retainQ;
ray_crosses_pd      =   Murat.PD.raysPeakDelay;
ray_crosses_Qc      =   Murat.Qc.raysQc;
ray_crosses_Q       =   Murat.Q.raysQ;
tCoda               =   Murat.Qc.tCoda;
outputSolverQc      =   struct;
outputSolverQ       =   struct;
exitFlagSolverQc    =   struct;
exitFlagSolverQ     =   struct;

% Paths and outputs
FPath               =   './';
outDirTXT           =   fullfile(FPath, FLabel, 'TXT');

lMF                 =   size(ray_crosses_pd);
nFreq               =   lMF(2);

% Preallocate
modv_pd             =   zeros(nxyzc,5,nFreq);
modv_Qc             =   zeros(nxyzc,10,nFreq);
modv_Q              =   zeros(nxyzc,10,nFreq);
const_Qc            =   zeros(size(rapsp));
residualQ           =   zeros(1,nFreq);
residualQc          =   zeros(1,nFreq);
dValueQc            =   zeros(1,nFreq);
dValueQ             =   zeros(1,nFreq);

% Diffusion constant - Wu 1985
D                   =   vS./Le_1/3./(1-B0);

% Loops over all frequencies and parameter models
for k = 1:nFreq
    
    D_k             =   D(k);
    if ~isempty(dampValueQc)
        dValueQc_k  =   dampValueQc(k);
    else
        dValueQc_k  =   [];
    end
    
    if ~isempty(dampValueQ)
        dValueQ_k  =   dampValueQ(k);
    else
        dValueQ_k  =   [];
    end
    
    if ~isempty(smoothValueQc)
        sValueQc_k  =   smoothValueQc(k);
    else
        sValueQc_k  =   [];
    end
    
    if ~isempty(smoothValueQ)
        sValueQ_k  =   smoothValueQ(k);
    else
        sValueQ_k  =   [];
    end
    
    % copy coordinates block
    modv_pd(:,1:3,k)=   modvP(:,1:3);
    modv_Qc(:,1:3,k)=   modvP(:,1:3);
    modv_Q(:,1:3,k) =   modvP(:,1:3);
    cf_k            =   cf(k);
    fld             =   sprintf('Hz%g', cf_k);
    fld             =   strrep(fld, '.', '_');
    fcName          =   strrep(num2str(cf_k),'.','_');

    % --- Peak delay standard regionalization ---
    rcpd_k          =   ray_crosses_pd(:,k);
    rtpd_k          =   retain_pd(:,k);
    Apd_k           =	Apd_i(rtpd_k,rcpd_k);
    lpdelta_k       =   lpdelta(rtpd_k,k);
    
    A_boxes         =   Apd_k>0;
    counts_box      =   sum(A_boxes,1);
    mpd             =   sum(A_boxes.*lpdelta_k,1)'./counts_box';

    mpd(isnan(mpd))     =   mean(mpd,'omitnan');
    modv_pd(rcpd_k,4,k) =   mpd;
    modv_pd(rcpd_k,5,k) =   counts_box;

    % --- Qc inversion ---
    rcQc_k          =   ray_crosses_Qc(:,k);
    rtQc_k          =   retain_Qc(:,k);
    Ac_k            =   Ac_i(rtQc_k,rcQc_k,k);
    RZZ_k           =   RZZ(rtQc_k,k);
    Qm_k            =   Qm(rtQc_k,k);
    coordPriorQc    =   modv_Qc(rcQc_k,1:3,k);
    FName           =   ['L-curve_Qc_' fcName '_Hz'];
    outDirFigure    =   fullfile(FPath,FLabel,'Tests/LCurve',FName);
    
    [sol,fval,exitflag,output,dValueQc_k]= Murat_inversionQc(Ac_k,Qm_k,...
        inversionMethod,iter,iterStall,RZZ_k,coordPriorQc,dValueQc_k,...
        sValueQc_k,plotI,outDirFigure);
    
    modv_Qc(rcQc_k,4,k)     =   sol.Qc;
    residualQc(k)           =   fval;
    fprintf(['Qc relative misfit reduction at ' num2str(cf_k)...
        'Hz is: %.4f\n'], fval);
    exitFlagSolverQc.(fld)  =   exitflag;
    outputSolverQc.(fld)    =   output;
    
    % --- Q inversion ---
    rcQ_k               =   ray_crosses_Q(:,k);
    rtQ_k               =   retain_Q(:,k);
    keepMask            =   rtQ_k & rtQc_k;
    A_k                 =   A_i(keepMask,rcQ_k);
    Qc_k                =   Qm(keepMask,k);
    luntot_k            =   luntot(keepMask);
    l                   =   luntot_k/1000;
    time0_k             =   time0(keepMask);
    rapsp_k             =   rapsp(keepMask,k);
    tCm                 =   tCoda(keepMask,k);
    coordPriorQ         =   modv_Q(rcQ_k,1:3,k);
    
    [d0, Q0,const_Qc_k] =   Murat_lsqlinQmean(tCm,tWm,Qc_k,cf_k,D_k,l,...
        time0_k,rapsp_k);
    
    FName               =   ['L-curve_Q_' fcName '_Hz'];
    outDirFigure        =   fullfile(FPath,FLabel,'Tests/LCurve',FName);
    
    [sol,fval,exitflag,output,dValueQ_k]  =   Murat_inversionQ(A_k,d0,...
        rapsp_k,inversionMethod,Q0,iter,iterStall,coordPriorQ,dValueQ_k,...
        sValueQ_k,plotI,outDirFigure);

    modv_Q(rcQ_k,4,k)       =   sol.Q;
    residualQ(k)            =   fval;
    fprintf(['Q relative misfit reduction at ' num2str(cf_k)...
        'Hz is: %.4f\n'], fval);
    exitFlagSolverQ.(fld)   =   exitflag;   
    outputSolverQ.(fld)     =   output;
    

    % --- Checkerboards and spike inputs and checkerboard inversion ---
    % --- Qc ---
    siz                     =   [nxc nyc nzc];
    I                       =   checkerBoard3D(siz,sizea);
    [checkInput,spikeInput] =  Murat_inputTesting(I,spike_o,spike_e,x,y,z);

    modv_Qc(checkInput==1,6,k)  =   min(Qm_k);
    modv_Qc(checkInput==0,6,k)  =   max(Qm_k);
    modv_Qc(:,8,k)              =   mean(Qm_k);
    modv_Qc(spikeInput,8,k)     =   spike_v;
    
    Qc_ch           =   modv_Qc(rcQc_k,6,k);
    
    re_checkQc      =   Ac_k*Qc_ch;
    [sol,~,~,~,~]   =   Murat_inversionQc(Ac_k,re_checkQc,...
        inversionMethod,iter,iterStall,RZZ_k,coordPriorQc,...
        dValueQc_k,sValueQc_k,plotI,outDirFigure);
    modv_Qc(rcQc_k,7,k) =   sol.Qc;
    
    % --- Q: reuse Qc checkerboard parameters ---
    modv_Q(:,6:8,k) =   modv_Qc(:,6:8,k);
    Q_ch            =   modv_Q(rcQ_k,6,k);
    re_checkQ       =   A_k*Q_ch;
    
    %Synthetic energy ratios
    synthEratio     =   exp(const_Qc_k + 2*pi*cf_k*re_checkQ)*2*pi*cf_k...
        ./l.^2;

    [d0c, Q0c, ~]   =   Murat_lsqlinQmean(tCm,tWm,Qc_k,cf_k,D_k,l,...
        time0_k,synthEratio);
    
    [sol,~,~,~,~]   =   Murat_inversionQ(A_k,d0c,rapsp_k,...
        inversionMethod,Q0c,iter,iterStall,coordPriorQ,dValueQ_k,...
        sValueQ_k,plotI,outDirFigure);
        
    modv_Q(rcQ_k,7,k)   =   sol.Q;

    % --- Spike inversion if requested ---
    if ~isempty(spike_o)
        Qc_sp           =   modv_Qc(rcQc_k,8,k);
        re_spikeQc      =   Ac_k*Qc_sp;
        [sol,~,~,~,~]   =   Murat_inversionQc(Ac_k,re_spikeQc,...
            inversionMethod,iter,iterStall,RZZ_k,coordPriorQc,...
        dValueQc_k,sValueQc_k,plotI,outDirFigure);
        modv_Qc(rcQc_k,9,k) =   sol.Qc;
        
        Q_sp            =   modv_Q(rcQ_k,8,k);
        re_spikeQ       =   A_k*Q_sp;
        
        %Synthetic energy ratios
        synthEratio     =   exp(const_Qc_k + 2*pi*cf_k*re_spikeQ)*...
            2*pi*cf_k./l.^2;
        
        [d0s, Q0s, ~]   =   Murat_lsqlinQmean(tCm,tWm,Qc_k,cf_k,D_k,l,...
            time0_k,synthEratio);
        
        [sol,~,~,~,~]   =   Murat_inversionQ(A_k,d0s,rapsp_k,...
            inversionMethod,Q0s,iter,iterStall,coordPriorQ,dValueQ_k,...
            sValueQ_k,plotI,outDirFigure);

        modv_Q(rcQ_k,9,k)   =   sol.Q;

    end

    % --- Save outputs for this frequency ---
    modv_pd_dd          =   modv_pd(:,:,k);
    modv_pd_dd(:,1:3)   =   DDcoordinates;
    modv_Qc_dd          =   modv_Qc(:,:,k);
    modv_Qc_dd(:,1:3)   =   DDcoordinates;
    modv_Q_dd           =   modv_Q(:,:,k);
    modv_Q_dd(:,1:3)    =   DDcoordinates;

    dValueQc(k)         =   dValueQc_k;
    dValueQ(k)          =   dValueQ_k;
    writematrix(modv_pd_dd,...
        fullfile(outDirTXT, ['peakdelay_' fcName '_Degrees_Hz.txt']));
    writematrix(modv_Qc_dd,...
        fullfile(outDirTXT, ['Qc_' fcName '_Degrees_Hz.txt']));
    writematrix(modv_Q_dd,...
        fullfile(outDirTXT, ['Q_' fcName '_Degrees_Hz.txt']));
end

% Save back to Murat
Murat.Qc.residualRedQc  =   residualQc;
Murat.Qc.dampingValueQc =   dValueQc;

Murat.Q.const_Qc        =   const_Qc;
Murat.Q.residualRedQ    =   residualQ;
Murat.Q.dampingValueQ   =   dValueQ;

Murat.PD.modvPeakDelay  =   modv_pd;
Murat.Qc.modvQc         =   modv_Qc;
Murat.Q.modvQ           =   modv_Q;

Murat.Qc.outputSolverQc =   outputSolverQc;
Murat.Q.outputSolverQ   =   outputSolverQ;

Murat.Qc.exitFlagSolverQc   =   exitFlagSolverQc;
Murat.Q.exitFlagSolverQ =   exitFlagSolverQ;
Murat.input.DiffConstant=   D;

writetable(muratHeader,fullfile(outDirTXT, 'DataHeaders.xlsx'));