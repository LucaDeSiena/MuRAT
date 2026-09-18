function Murat      =   Murat_dataParallelized(Murat)
% MEASURES Qc, peak-delay and Q for each seismic trace located in a folder.

% Refactored input unpacking and preallocation
s                   =   Murat.input;
FLabel              =   s.label;
workers             =   s.workers; %#ok<NASGU>

%  Basic lists and sizes
listSac             =   s.listSac;
nData               =   numel(listSac);
sacHeader           =   s.sacHeader;

components          =   s.components;
sections            =   s.sections;
sections(3)         =   sections(3)/1000;
modv                =   s.modv;
nModelParams        =   size(modv(:,1),1);

% Grid step and model grid coordinates
gridStepZ           =   (modv(2,3) - modv(1,3))/2;
gridStep            =   [s.gridStepX/2, s.gridStepY/2, gridStepZ];
modvQc              =   [modv(:,1) + gridStep(1)...
    modv(:,2)+gridStep(2) modv(:,3)+gridStep(3)];

gridD               =   s.gridPropagation;
pvel                =   s.pvel;

% Frequencies
cf                  =   s.centralFrequency;
nCF                 =   numel(cf);

% Event/origin and timing
origin              =   s.origin;
originTime          =   s.originTime;
PTime               =   s.PTime;
STime               =   s.STime;
PorS                =   s.POrS;
startLapse          =   s.startLapseTime;

% Velocities and travel times
vP                  =   s.averageVelocityP;
vS                  =   s.averageVelocityS;
maxtpde             =   s.maximumPeakDelay;
maxTravel           =   s.maxtravel;
minTravel           =   s.mintravel;
lapseTimeMethod     =   s.lapseTimeMethod;

% Coda parameters
tCodaWindow         =   s.codaWindow;
kernelThresh        =   s.kernelTreshold;
B0                  =   s.albedo;
Le                  =   s.iExtinctionLength;
bodyWindow          =   s.bodyWindow;
startNoise          =   s.startNoise;
QcMeasure           =   s.QcMeasurement;


% Preallocate outputs
locationDeg         =   zeros(nData,6); 
locationM           =   zeros(nData,6); 
theoreticalTime     =   zeros(nData,1); 
totalLengthRay      =   zeros(nData,1);
peakDelay           =   nan(nData,nCF);
inverseQc           =   nan(nData,nCF);
uncertaintyQc       =   nan(nData,nCF);
energyRatioBCo      =   nan(nData,nCF);
energyRatioCNo      =   nan(nData,nCF);
raysPlot            =   zeros(100,5,nData);
tCoda               =   zeros(nData,nCF);

inversionMatrixPD   =   zeros(nData,nModelParams);
inversionMatrixQ    =   zeros(nData,nModelParams);
inversionMatrixQc   =   zeros(nData,nModelParams,nCF);
rayCrossing         =   zeros(nData,nModelParams);

storeFolderK        =   'Kernels';
       
% Count waveforms that must be eliminated because of peak-delay contraints
% on peak delays and coda attenuation.
countTrash          =   0;
%=========================================================================
SAChdrList = cell(nData,1);
for k = 1:nData
    fdl = sprintf('SAC_%g', k);
    SAChdrList{k} = sacHeader.(fdl);
end

parfor (i = 1:nData,workers)
    
    % Progress every 100 traces
    if isequal(mod(i,100),0)

        fprintf('Waveform number is %d\n', i);
        
    end
    
    listSac_i       =   listSac{i};
    SAChdr_i        =   SAChdrList{i};
    srate_i         =   1/SAChdr_i.times.delta;
    
    % Envelope and spectral preprocessing
    [tempis,sp_i]   =   Murat_envelope(cf,listSac_i);
    
    % Locations (deg and meters)
    [locationDeg_i, locationM_i]    =   Murat_location(origin,SAChdr_i);
    locationDeg(i,:)                =   locationDeg_i;
    
    % Checks direct-wave picking 
    [cursorPick_i, pktime_i, v_i]   =   Murat_picking(tempis,...
        PTime,STime,PorS,vP,vS,srate_i,listSac_i,SAChdr_i);

    % Missing origin-time corrections
    [theoreticalTime_i, originTime_i]   =...
        Murat_originTime(pktime_i,originTime,v_i,locationM_i,SAChdr_i,i);

    % Peak-delay search window and value
    cursorPeakDelay_i           =...
        Murat_peakDelayCheck(tempis,cursorPick_i,maxtpde,srate_i);
    
    peakDelay_i                 =...
        Murat_peakDelay(sp_i,cursorPick_i,srate_i,cursorPeakDelay_i);
    
    % Decide if ray-dependent computations are required for this component
    doRays                      =   recognizeComponents(i, components);
    
    if doRays
        % All the ray-dependent parameters   
        [Apd_i, AQ_i, totalLengthRay_i, raysPlot_i, rayCrossing_i]=...                                        
            Murat_rays(modv,gridD,pvel,locationM_i);        
        inversionMatrixPD(i,:)  =   Apd_i;
        inversionMatrixQ(i,:)   =   AQ_i;
        totalLengthRay(i,1)     =   totalLengthRay_i;
        raysPlot(:,:,i)         =   raysPlot_i;
        rayCrossing(i,:)        =   rayCrossing_i;
    end     
    % 
    % % Diffusion equation
    % [storeD,residualD,outputSolverD,exitFlagSolverD] =...
    %         Murat_diffusion(cf,sp_i,startLapse,tCodaWindow,totalLengthRay_i);
    
    % Coda window, checks, and lapse time
    [tCoda_i, cursorCodaStart_i,cursorCodaEnd_i] = Murat_codaCheck(...
        originTime_i,pktime_i,srate_i,startLapse,tCodaWindow,tempis,...
        peakDelay_i,lapseTimeMethod);
    
    % Reject traces that don't meet coda/travel requirements
    codaSamples                 =  cursorCodaEnd_i - cursorCodaStart_i;
    if mean(codaSamples) < (tCodaWindow*srate_i) - 2 || ...
            (pktime_i-originTime_i) > maxTravel || ...
            (pktime_i-originTime_i) < minTravel

        locationM(i,:)          =   locationM_i;
        theoreticalTime(i,1)    =   theoreticalTime_i;
        tCoda(i,:)              =   tCoda_i;
        countTrash              =   countTrash +1;
        continue
    end
    
    % Measures Qc and its uncertainty
    [inverseQc_i, uncertaintyQc_i]  =   Murat_Qc(cf,...
        sp_i,cursorCodaStart_i,cursorCodaEnd_i,tCoda_i,srate_i,QcMeasure);
    
    % Decide if kernel/inversion matrix for Qc is needed
    doKernels                   =   doRays;
    
    if doKernels
        % Compute and plots kernels and coda matrix
        TCodaStart              =   tCoda_i;
        TCodaEnd                =   tCoda_i+tCodaWindow;
        [KgridS, rgridS]        =...
            Murat_kernels(TCodaStart,locationM_i(1:3),locationM_i(4:6),...
            modvQc,vS,kernelThresh,B0,Le,lapseTimeMethod);
        [KgridE, rgridE]        =...
            Murat_kernels(TCodaEnd,locationM_i(1:3),locationM_i(4:6),...
            modvQc,vS,kernelThresh,B0,Le,lapseTimeMethod);
        
        for j = 1:nCF

            cf_k    =   cf(j);
            fcName  =   num2str(cf_k);
            if find(fcName == '.')
                fcName(fcName == '.') =   '_';
            end
            FName   =   ['Kernel_' fcName '_Hz'];
            KS_i    =   KgridS(:,:,j);
            KE_i    =   KgridE(:,:,j);
            p       =   fullfile('./', FLabel, storeFolderK, FName);
        
            AQc_i   =   Murat_codaMatrix(modvQc,KS_i,rgridS,KE_i,...
                rgridE,i,origin,sections,FName,p);            
        
            inversionMatrixQc(i,:,j)  =   AQc_i;

        end
    end
                
    % Measure body/coda energy ratios
    [energyRatioBodyCoda_i,energyRatioCodaNoise_i]=...
        Murat_body(bodyWindow,startNoise,srate_i,sp_i,cursorPick_i,...
        cursorCodaStart_i,cursorCodaEnd_i);
    
    % Save results for this trace
    locationM(i,:)              =   locationM_i;
    theoreticalTime(i,1)        =   theoreticalTime_i;
    peakDelay(i,:)              =   peakDelay_i;
    inverseQc(i,:)              =   inverseQc_i; 
    uncertaintyQc(i,:)          =   uncertaintyQc_i; 
    energyRatioBCo(i,:)         =   energyRatioBodyCoda_i; 
    energyRatioCNo(i,:)         =   energyRatioCodaNoise_i;
    tCoda(i,:)                  =   tCoda_i;
    
end

% Assign to Murat structure
Murat.rays.locationsDeg         =   locationDeg;
Murat.rays.locationsM           =   locationM;
Murat.rays.theoreticalTime      =   theoreticalTime;
Murat.PD.peakDelay              =   peakDelay;
Murat.PD.inversionMatrixPeakDelay   =   inversionMatrixPD;
Murat.Q.inversionMatrixQ        =   inversionMatrixQ;
Murat.rays.totalLengthRay       =   totalLengthRay;
Murat.rays.raysPlot             =   raysPlot;
Murat.rays.rayCrossing          =   sum(rayCrossing);
Murat.Qc.inverseQc              =   inverseQc; 
Murat.Qc.uncertaintyQc          =   uncertaintyQc; 
Murat.Qc.inversionMatrixQc      =   inversionMatrixQc;
Murat.Q.energyRatioBodyCoda     =   energyRatioBCo; 
Murat.Q.energyRatioCodaNoise    =   energyRatioCNo;
Murat.Qc.tCoda                  =   tCoda;

Murat                           =   Murat_selection(Murat);
ratio                           =   countTrash/nData*(100);
fprintf('Ratio of recordings removed due to Qc or PD inputs: %.2f%%\n',...
    ratio)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doR                   =  recognizeComponents(index,components)
% LOGICAL to decide if forward model is necessary depending in waveform
%   number (index) and number of components.

if components == 1
    doR = true;
elseif components == 2 || components == 3
    % Only process every components-th trace starting at 1
    doR = mod(index, components) == 1;
else
    doR = false;
end
end
