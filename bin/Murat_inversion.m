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
rapE                =   Murat.Q.energyRatioBodyCoda;
retainPD            =   Murat.PD.retainPeakDelay;
retainQc            =   Murat.Qc.retainQc;
retainQ             =   Murat.Q.retainQ;
rayCrossesPD        =   Murat.PD.raysPeakDelay;
rayCrossesQc        =   Murat.Qc.raysQc;
rayCrossesQ         =   Murat.Q.raysQ;
tCoda               =   Murat.Qc.tCoda;
outputSolverQc      =   struct;
outputSolverQ       =   struct;
exitFlagSolverQc    =   struct;
exitFlagSolverQ     =   struct;

% Paths and outputs
FPath               =   './';
outDirTXT           =   fullfile(FPath, FLabel, 'TXT');
lMF                 =   size(rayCrossesPD);
nFreq               =   lMF(2);

% Preallocate
modvPD              =   zeros(nxyzc,5,nFreq);
modvQc              =   zeros(nxyzc,9,nFreq);
modvQ               =   zeros(nxyzc,9,nFreq);
residualQ           =   zeros(1,nFreq);
residualQc          =   zeros(1,nFreq);
dValueQc            =   zeros(1,nFreq);
dValueQ             =   zeros(1,nFreq);
cConstants          =   zeros(2,nFreq);

% Diffusion constant - Wu 1985, Wu and Aki, 1988.
D                   =   vS/3./Le_1./B0;

fc_names    = cell(1,nFreq);
fld_names   = cell(1,nFreq);
outDirFigureBase    = fullfile(FPath,FLabel,'Tests','LCurve');
for k = 1:nFreq
    cf_k    = cf(k);
    fc_names{k}     = strrep(num2str(cf_k),'.','_');
    tmp     = sprintf('Hz%g', cf_k);
    fld_names{k}    = strrep(tmp,'.','_');
end

% reuse buffers to avoid reallocations
siz     = [nxc nyc nzc];
I   =   []; checkInput    =     []; spikeInput    =     [];

% Precompute coordinate block once and reuse
coords_block = modvP(:,1:3);

% Loops over all frequencies and parameter models
for k = 1:nFreq
    
    D_k             =   D(k);
    if ~isempty(dampValueQc), dValueQc_k  =   dampValueQc(k);
    else, dValueQc_k    =   [];   end
    
    if ~isempty(dampValueQ), dValueQ_k    =   dampValueQ(k);
    else, dValueQ_k     =   [];   end
    
    if ~isempty(smoothValueQc), sValueQc_k=   smoothValueQc(k);
    else, sValueQc_k    =   []; end
    
    if ~isempty(smoothValueQ), sValueQ_k  =   smoothValueQ(k);
    else, sValueQ_k     =   []; end
    
    modvPD(:,1:3,k) =   coords_block;
    modvQc(:,1:3,k) =   coords_block;
    modvQ(:,1:3,k)  =   coords_block;
    cf_k            =   cf(k);
    fcName          =   fc_names{k};
    fld = fld_names{k};
    outDirFigure = fullfile(outDirFigureBase,['L-curve_Qc_' fcName '_Hz']);

    % --- Peak delay standard regionalization ---
    rcpd_k          =   rayCrossesPD(:,k);
    rtpd_k          =   retainPD(:,k);
    Apd_k           =	Apd_i(rtpd_k,rcpd_k);
    lpdelta_k       =   lpdelta(rtpd_k,k);
    
    A_boxes         =   Apd_k>0;
    counts_box      =   sum(A_boxes,1);
    mpd             =   sum(A_boxes.*lpdelta_k,1)'./counts_box';

    mpd(isnan(mpd))     =   mean(mpd,'omitnan');
    modvPD(rcpd_k,4,k)  =   mpd;
    modvPD(rcpd_k,5,k)  =   counts_box;

    % --- Qc inversion ---
    rcQc_k          =   rayCrossesQc(:,k);
    rtQc_k          =   retainQc(:,k);
    Ac_k            =   Ac_i(rtQc_k,rcQc_k,k);
    RZZ_k           =   RZZ(rtQc_k,k);
    Qm_k            =   Qm(rtQc_k,k);
    coordPriorQc    =   modvQc(rcQc_k,1:3,k);
    
    [sol,fval,exitflag,output,dValueQc_k]= Murat_inversionQc(Ac_k,Qm_k,...
        inversionMethod,iter,iterStall,RZZ_k,coordPriorQc,dValueQc_k,...
        sValueQc_k,plotI,outDirFigure);
    
    modvQc(rcQc_k,4,k)      =   sol.Qc;
    residualQc(k)           =   fval;
    fprintf(['Qc relative misfit reduction at ' num2str(cf_k)...
        'Hz is: %.4f\n'], fval);
    exitFlagSolverQc.(fld)  =   exitflag;
    outputSolverQc.(fld)    =   output;
    
    % --- Q inversion ---
    rcQ_k               =   rayCrossesQ(:,k);
    rtQ_k               =   retainQ(:,k);
    keepMask            =   rtQ_k & rtQc_k;
    A_k                 =   A_i(keepMask,rcQ_k);
    Qc_k                =   Qm(keepMask,k);
    luntot_k            =   luntot(keepMask);
    l                   =   luntot_k/1000;
    time0_k             =   time0(keepMask);
    rapE_k              =   rapE(keepMask,k);
    tCm                 =   tCoda(keepMask,k);
    coordPriorQ         =   modvQ(rcQ_k,1:3,k);
    te                  =   tCm+tWm;

    cCoda               =   1.5/2/pi/cf_k*log(4*pi*D_k)-...
                            1/2/pi/cf_k*...
                            log((te).^(-1.5).*...
                            exp(-l.^2/4/D_k./te...
                            -2*pi*cf_k*Qc_k.*te)...
                            -tCm.^(-1.5).*...
                            exp(-l.^2/4/D_k./tCm...
                            -2*pi*cf_k*Qc_k.*tCm));

    [d0, Q0, mDest]     =   Murat_lsqlinQmean(cf_k,l,time0_k,rapE_k,...
        Qc_k, tCm, te);

    cEst                =   Q0(1,1);

    cConstants(1,k)     =   mean(real(cCoda));
    cConstants(2,k)     =   cEst;
    
    fprintf('A priori and estimated coda constants are %.4g and %.4g at %g Hz\n',...
        mean(real(cCoda)), cEst, cf_k);
    
    fprintf('A priori and estimated diffusion constants are %.4g and %.4g at %g Hz\n',...
        D_k, mDest, cf_k);

    averageQcodaC       =   Q0;
    averageD            =   mDest;
    outDirFigure = fullfile(outDirFigureBase, ['L-curve_Q_' fcName '_Hz']);
    
    [sol,fval,exitflag,output,dValueQ_k]  =   Murat_inversionQ(A_k,d0,...
        rapE_k,inversionMethod,Q0,iter,iterStall,coordPriorQ,dValueQ_k,...
        sValueQ_k,plotI,outDirFigure);

    modvQ(rcQ_k,4,k)        =   sol.Q;
    residualQ(k)            =   fval;
    fprintf(['Q relative misfit reduction at ' num2str(cf_k)...
        'Hz is: %.4f\n'], fval);
    exitFlagSolverQ.(fld)   =   exitflag;   
    outputSolverQ.(fld)     =   output;

    % --- Checkerboards and spike inputs and checkerboard inversion ---
    % --- Qc ---
    if isempty(I)
        I = checkerBoard3D(siz, sizea);
        [checkInput, spikeInput] =   Murat_inputTesting(I, spike_o,...
            spike_e, x, y, z);
    end
    
    Qm_k_all = Qm(rtQc_k,k);  % local
    
    modvQc(checkInput==1,6,k)  =   min(Qm_k_all);
    modvQc(checkInput==0,6,k)  =   max(Qm_k_all);
    modvQc(:,8,k)              =   mean(Qm_k_all);
    if ~isempty(spikeInput)
        modvQc(spikeInput,8,k) = spike_v;
    end
    
    Qc_ch           =   modvQc(rcQc_k,6,k);
    re_checkQc      =   Ac_k*Qc_ch;
    [sol,~,~,~,~]   =   Murat_inversionQc(Ac_k,re_checkQc,...
        inversionMethod,iter,iterStall,RZZ_k,coordPriorQc,...
        dValueQc_k,sValueQc_k,plotI,outDirFigure);
    modvQc(rcQc_k,7,k)  =   sol.Qc;
    
    % --- Q: reuse Qc checkerboard parameters ---
    modvQ(:,6:8,k)  =   modvQc(:,6:8,k);
    Q_ch            =   modvQ(rcQ_k,6,k);
    re_checkQ       =   A_k*Q_ch;
    
    % Synthetic energy ratios
    synthEratio     =   exp(cEst + 2*pi*cf_k*re_checkQ)*2*pi*cf_k...
        ./l.^2;

    [d0c, Q0c, ~]   =   Murat_lsqlinQmean(cf_k,l,time0_k,synthEratio,...
        Qc_k, tCm, te);
    
    [sol,~,~,~,~]   =   Murat_inversionQ(A_k,d0c,rapE_k,...
        inversionMethod,Q0c,iter,iterStall,coordPriorQ,dValueQ_k,...
        sValueQ_k,plotI,outDirFigure);
        
    modvQ(rcQ_k,7,k)=   sol.Q;

    % --- Spike inversion if requested ---
    if ~isempty(spike_o)
        Qc_sp           =   modvQc(rcQc_k,8,k);
        re_spikeQc      =   Ac_k*Qc_sp;
        [sol,~,~,~,~]   =   Murat_inversionQc(Ac_k,re_spikeQc,...
            inversionMethod,iter,iterStall,RZZ_k,coordPriorQc,...
        dValueQc_k,sValueQc_k,plotI,outDirFigure);
        modvQc(rcQc_k,9,k) =   sol.Qc;
        
        Q_sp            =   modvQ(rcQ_k,8,k);
        re_spikeQ       =   A_k*Q_sp;
        
        %Synthetic energy ratios
        synthEratio     =   exp(cEst + 2*pi*cf_k*re_spikeQ)*...
            2*pi*cf_k./l.^2;
        
        [d0s, Q0s, ~]   =  Murat_lsqlinQmean(cf_k,l,time0_k,synthEratio,...
            Qc_k, tCm, te);
        
        [sol,~,~,~,~]   =   Murat_inversionQ(A_k,d0s,rapE_k,...
            inversionMethod,Q0s,iter,iterStall,coordPriorQ,dValueQ_k,...
            sValueQ_k,plotI,outDirFigure);

        modvQ(rcQ_k,9,k)=   sol.Q;

    end

    % --- Save outputs for this frequency ---
    modv_pd_dd          =   modvPD(:,:,k);
    modv_pd_dd(:,1:3)   =   DDcoordinates;
    modv_Qc_dd          =   modvQc(:,:,k);
    modv_Qc_dd(:,1:3)   =   DDcoordinates;
    modv_Q_dd           =   modvQ(:,:,k);
    modv_Q_dd(:,1:3)    =   DDcoordinates;
    dValueQc(k)         =   dValueQc_k;
    dValueQ(k)          =   dValueQ_k;

    headerPD            =   {'Lat','Lon','Depth','ΔlogPD', 'HitC'};
    C                   =   [headerPD; num2cell(modv_pd_dd)];
    fname               =...
        fullfile(outDirTXT, ['peakdelay_' fcName '_Hz.txt']);
    writecell(C, fname, 'Delimiter', '\t', 'Encoding', 'UTF-8');
    
    headerQc            =...
        {'Lat','Lon','Depth','Qc', 'HitC','Check In','Check Out',...
        'Spike In','Spike Out'};
  C                   =   [headerQc; num2cell(modv_Qc_dd)];
    fname               =...
        fullfile(outDirTXT, ['Qc_' fcName '_Hz.txt']);
    writecell(C, fname, 'Delimiter', '\t', 'Encoding', 'UTF-8');
   
    
    headerQ            =...
        {'Lat','Lon','Depth','Q', 'HitC','Check In','Check Out',...
        'Spike In','Spike Out'};
    C                   =   [headerQ; num2cell(modv_Q_dd)];
    fname               =...
        fullfile(outDirTXT, ['Q_' fcName '_Hz.txt']);
    writecell(C, fname, 'Delimiter', '\t', 'Encoding', 'UTF-8');
    
end

% Save back to Murat
Murat.input.DiffConstant    =   D;

Murat.PD.modvPeakDelay      =   modvPD;

Murat.Qc.modvQc             =   modvQc;
Murat.Qc.residualRedQc      =   residualQc;
Murat.Qc.dampingValueQc     =   dValueQc;
Murat.Qc.outputSolverQc     =   outputSolverQc;
Murat.Qc.exitFlagSolverQc   =   exitFlagSolverQc;

Murat.Q.residualRedQ        =   residualQ;
Murat.Q.dampingValueQ       =   dValueQ;
Murat.Q.modvQ               =   modvQ;
Murat.Q.outputSolverQ       =   outputSolverQ;
Murat.Q.exitFlagSolverQ     =   exitFlagSolverQ;
Murat.Q.codaConstants       =   cConstants;
Murat.Q.averageQ            =   averageQcodaC;
Murat.Q.estimatedDiffusion  =   averageD;
end