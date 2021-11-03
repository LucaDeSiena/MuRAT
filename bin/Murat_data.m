function Murat                          =   Murat_data(Murat)
% MEASURES Qc, peak-delay and Q for each seismic trace located in a folder.
%   This code is a collection of functions that do all the necessary.

% Inputs
listSac                                 =   Murat.input.listSac;
lengthData                              =   length(listSac);

compon                                  =   Murat.input.components;

modv                                    =   Murat.input.modv;
lengthParameterModel                    =   length(modv(:,1));
gridStep                                =   [Murat.input.gridStepX/2 ...
    Murat.input.gridStepY/2 (modv(2,3) - modv(1,3))/2];
modvQc                                  =   [modv(:,1) + gridStep(1)...
    modv(:,2)+gridStep(2) modv(:,3)+gridStep(3)];

gridD                                   =   Murat.input.gridPropagation;
pvel                                    =   Murat.input.pvel;

cf                                      =   Murat.input.centralFrequency;
lcf                                     =   length(cf);

origin                                  =   Murat.input.origin;
originTime                              =   Murat.input.originTime;
PTime                                   =   Murat.input.PTime;
STime                                   =   Murat.input.STime;
PorS                                    =   Murat.input.POrS;
tCm                                     =   Murat.input.startLapseTime;
vP                                      =   Murat.input.averageVelocityP;
vS                                      =   Murat.input.averageVelocityS;
maxtpde                                 =   Murat.input.maximumPeakDelay;
tWm                                     =   Murat.input.codaWindow;
sped                                    =   Murat.input.spectralDecay;
kT                                      =   Murat.input.kernelTreshold;
B0                                      =   Murat.input.albedo;
Le1                                     =   Murat.input.extinctionLength;
bodyWindow                              =   Murat.input.bodyWindow;
startNoise                              =   Murat.input.startNoise;
QcM                                     =   Murat.input.QcMeasurement;
lapseTimeMethod                         =   Murat.input.lapseTimeMethod;
maxtravel                               =   Murat.input.maxtravel;

% Set up variables to save
locationDeg                             =   zeros(lengthData,6); 
locationM                               =   zeros(lengthData,6); 
theoreticalTime                         =   zeros(lengthData,1); 
totalLengthRay                          =   zeros(lengthData,1);
peakDelay                               =   zeros(lengthData,lcf); 
inverseQc                               =   zeros(lengthData,lcf); 
uncertaintyQc                           =   zeros(lengthData,lcf); 
energyRatioBodyCoda                     =   zeros(lengthData,lcf); 
energyRatioCodaNoise                    =   zeros(lengthData,lcf);
raysPlot                                =   zeros(100,5,lengthData);
tCoda                                   =   zeros(lengthData,lcf);

inversionMatrixPeakDelay                =...
    zeros(lengthData,lengthParameterModel);
inversionMatrixQ                        =...
    zeros(lengthData,lengthParameterModel);
inversionMatrixQc                       =...
    zeros(lengthData,lengthParameterModel);
rayCrossing                             =...
    zeros(lengthData,lengthParameterModel);
%=========================================================================
count_trash = 0;
for i = 1:lengthData
    
    if isequal(mod(i,1000),0)
        
        disp(['Waveform number is ', num2str(i)])
        
    end
    
    listSac_i                           =   listSac{i};
    
    % Calculates envelopes
    [tempis,sp_i,SAChdr_i,srate_i]      =...
        Murat_envelope(cf,listSac_i);
    
    % Set earthquake and stations locations in degrees or meters
    [locationDeg_i, locationM_i]        =...
        Murat_location(origin,SAChdr_i);
    locationDeg(i,:)                    =   locationDeg_i;
    
    % Checks direct-wave picking on the trace and outputs it 
    [cursorPick_i, pktime_i, v_i]  =...
        Murat_picking(tempis,PTime,STime,PorS,vP,vS,srate_i,listSac_i,...
        SAChdr_i);

    % Conditions in case the zero time is missing in the header
    [theoreticalTime_i, originTime_i]   =...
        Murat_originTime(pktime_i,originTime,v_i,locationM_i,SAChdr_i);

    % Calculates the window where to search for peak delay
    cursorPeakDelay_i                   =...
        Murat_peakDelayCheck(tempis,cursorPick_i,maxtpde,srate_i);

    % Calculates peak delay time
    peakDelay_i                         =...
        Murat_peakDelay(sp_i,cursorPick_i,srate_i,cursorPeakDelay_i);
    
    % Calculates rays for the right component    
    calculateRays                       =   recognizeComponents(i,compon);
    
    if calculateRays
        
        % All the ray-dependent parameters   
        [Apd_i, AQ_i, totalLengthRay_i, raysPlot_i, rayCrossing_i]...
                                        =...
            Murat_rays(modv,gridD,pvel,locationM_i);
        
        inversionMatrixPeakDelay(i,:)   =   Apd_i;
        inversionMatrixQ(i,:)           =   AQ_i;
        totalLengthRay(i,1)             =   totalLengthRay_i;
        raysPlot(:,:,i)                 =   raysPlot_i;
        rayCrossing(i,:)                =   rayCrossing_i;
    end
                
    % Sets the lapse time
    [tCoda_i, cursorCodaStart_i, cursorCodaEnd_i]=...
        Murat_codaCheck(originTime_i,pktime_i,srate_i,tCm,tWm,tempis,...
        peakDelay_i,lapseTimeMethod);
    
    if (cursorCodaEnd_i -cursorCodaStart_i)< (tWm*srate_i)-2 || ...
            (pktime_i-originTime_i)>maxtravel

        locationM(i,:)                      =   locationM_i;
        theoreticalTime(i,1)                =   theoreticalTime_i;
        peakDelay(i,:)                      =   NaN;
        inverseQc(i,:)                      =   NaN;
        uncertaintyQc(i,:)                  =   NaN;
        energyRatioBodyCoda(i,:)            =   NaN;
        energyRatioCodaNoise(i,:)           =   NaN;
        tCoda(i,:)                          =   tCoda_i;
        
        count_trash = count_trash +1;
        continue
    end
    
    % Measures Qc and its uncertainty
    [inverseQc_i, uncertaintyQc_i]      =   Murat_Qc(cf,sped,...
        sp_i,cursorCodaStart_i,cursorCodaEnd_i,tCoda_i,srate_i,QcM);
    
    % Decide if you calculate kernels
    calculateKernels                    =   recognizeComponents(i,compon);
    
    if calculateKernels
        
        % Calculates kernels
        [K_grid, r_grid]                =...
            Murat_kernels(tCoda_i+tWm/2,locationM_i(1:3),...
            locationM_i(4:6),modvQc,vS,kT,B0,Le1,lapseTimeMethod);
        
        % Calculates matrix
        AQc_i                           =...
            Murat_codaMatrix(modvQc,K_grid,r_grid,0,[],[]);
            
        inversionMatrixQc(i,:)          =   AQc_i;
        
    end
                
    % Measures Q
    [energyRatioBodyCoda_i,...
        energyRatioCodaNoise_i]         =   Murat_body(bodyWindow,...
        startNoise,srate_i,sp_i,cursorPick_i,...
        cursorCodaStart_i,cursorCodaEnd_i);
    
    % Saving
    locationM(i,:)                      =   locationM_i;
    theoreticalTime(i,1)                =   theoreticalTime_i;
    peakDelay(i,:)                      =   peakDelay_i;
    inverseQc(i,:)                      =   inverseQc_i; 
    uncertaintyQc(i,:)                  =   uncertaintyQc_i; 
    energyRatioBodyCoda(i,:)            =   energyRatioBodyCoda_i; 
    energyRatioCodaNoise(i,:)           =   energyRatioCodaNoise_i;
    tCoda(i,:)                          =   tCoda_i;
    
end

% Setting up the final data vectors and matrices with checks on values
Murat.data.locationsDeg                 =   locationDeg;
Murat.data.locationsM                   =   locationM;
Murat.data.theoreticalTime              =   theoreticalTime;
Murat.data.peakDelay                    =   peakDelay;
Murat.data.inversionMatrixPeakDelay     =   inversionMatrixPeakDelay;
Murat.data.inversionMatrixQ             =   inversionMatrixQ;
Murat.data.totalLengthRay               =   totalLengthRay;
Murat.data.raysPlot                     =   raysPlot;
Murat.data.rayCrossing                  =   sum(rayCrossing);
Murat.data.inverseQc                    =   inverseQc; 
Murat.data.uncertaintyQc                =   uncertaintyQc; 
Murat.data.inversionMatrixQc            =   inversionMatrixQc;
Murat.data.energyRatioBodyCoda          =   energyRatioBodyCoda; 
Murat.data.energyRatioCodaNoise         =   energyRatioCodaNoise;
Murat.data.tCoda                        =   tCoda;

Murat                                   =   Murat_selection(Murat);

ratio                                   =   count_trash/i*(100);
disp(['trash ', num2str(ratio)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calculateValue                 =...
    recognizeComponents(index,components)
% LOGICAL to decide if forward model is necessary depending in waveform
%   number (index) and number of components.

calculateValue                          =   isequal(components,1) ||...
    (isequal(components,2) || isequal(components,3)) &&...
    isequal(mod(index,components),1);
