%% Seismic attributes for peak delay and Qc imaging
function Murat                          =   Murat_dataParallelized(Murat)
%MEASURES Qc, peak-delay and Q for each seismic trace located in a folder.
% It only accepts SAC files.
% In the case of more than one component the rays must be ordered in the
% folder starting with the WE component, then SN and Z for each ray.
%
% The program uses the information in the SAC header, which MUST include
% peaking of the direct phase of interest (preferred is marker "a" for P
% and "t0" for S), as well as event and station locations.

listSac                                 =   Murat.input.listSac;
lengthData                              =   length(listSac);
compon                                  =   Murat.input.components;
lengthParameterModel                    =   length(Murat.input.modv(:,1));
cf                                      =   Murat.input.centralFrequency;
lcf                                     =   length(cf);

%Set up variables to save
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

inversionMatrixPeakDelay                =...
    zeros(lengthData,lengthParameterModel);
inversionMatrixQ                        =...
    zeros(lengthData,lengthParameterModel);
inversionMatrixQc                       =...
    zeros(lengthData,lengthParameterModel);
rayCrossing                             =...
    zeros(lengthData,lengthParameterModel);

%=========================================================================
parfor i = 1:lengthData %loop through source-station pairs
    
    if isequal(mod(i,1000),0)
        
        disp(['Waveform number is ', num2str(i)])
        
    end
    
    %% OPERATIONS ON WAVEFORM
    listSac_i                           =   listSac{i};
    [sp_i,SAChdr_i,srate_i]             =...
        Murat_envelope(Murat,listSac_i);
    
    % In case it has been precalculated from external files.
    [locationDeg_i, locationM_i]        =   Murat_location(Murat,SAChdr_i);
    locationDeg(i,:)                    =   locationDeg_i;
    
    [theoreticalTime_i, tCoda_i, cursorPick_i, cursorPeakDelay_i,...
        cursorCodaStart_i, cursorCodaEnd_i]...
                                        =...
        Murat_times(Murat,locationM_i,listSac_i);
    
    %% OPERATIONS TO MEASURE AND MODEL PEAK DELAYS + INVERSION MATRIX Q
    peakDelay_i                         =...
        Murat_peakDelay(sp_i,cursorPick_i,cursorPeakDelay_i,srate_i);
    
    calculateRays                       =   recognizeComponents(i,compon);
    
    if calculateRays
        
        [Apd_i, AQ_i, totalLengthRay_i, raysPlot_i,...
            rayCrossing_i]              =...
            Murat_rays(Murat,locationM_i);
        
        inversionMatrixPeakDelay(i,:)   =   Apd_i;
        inversionMatrixQ(i,:)           =   AQ_i;
        totalLengthRay(i,1)             =   totalLengthRay_i;
        raysPlot(:,:,i)                 =   raysPlot_i;
        rayCrossing(i,:)                =   rayCrossing_i;
      
    end
                
    %% OPERATIONS TO MEASURE AND MODEL Qc
    [inverseQc_i,uncertaintyQc_i]       =...
        Murat_Qc(Murat,sp_i,cursorCodaStart_i,...
        cursorCodaEnd_i,tCoda_i,srate_i);
    
    calculateKernels                    =...
        recognizeComponents(i,compon);
    
    if calculateKernels
        
            AQc_i                       =...
                Murat_codaMatrix(Murat,0,tCoda_i,locationM_i);
            
            inversionMatrixQc(i,:)      =   AQc_i;
    
    end
                
    %% OPERATIONS TO MEASURE  Q
    [energyRatioBodyCoda_i,...
        energyRatioCodaNoise_i]         =   Murat_body(Murat,...
        srate_i,sp_i,cursorPick_i,cursorCodaStart_i,cursorCodaEnd_i);
    
    %% SAVING
    locationM(i,:)                      =   locationM_i;
    theoreticalTime(i,1)                =   theoreticalTime_i;
    peakDelay(i,:)                      =   peakDelay_i;
    inverseQc(i,:)                      =   inverseQc_i; 
    uncertaintyQc(i,:)                  =   uncertaintyQc_i; 
    energyRatioBodyCoda(i,:)            =   energyRatioBodyCoda_i; 
    energyRatioCodaNoise(i,:)           =   energyRatioCodaNoise_i;
    
end

%% Setting up the final data vectors and matrices with checks on values
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

%% SELECTS DATA FOR INVERSION
Murat                                   =   Murat_selection(Murat);


function calculateValue                 =...
    recognizeComponents(index,components)
%LOGICAL to decide if forward model is necessary depending in waveform
%number (index) and number of components.

calculateValue                          =   isequal(components,1) ||...
    (isequal(components,2) || isequal(components,3)) &&...
    isequal(mod(index,components),1);
