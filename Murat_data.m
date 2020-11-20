%% Seismic attributes for peak delay and Qc imaging
function Murat                          =   Murat_data(Murat)
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
origin                                  =   Murat.input.origin;
lengthParameterModel                    =   length(Murat.input.modv(:,1));

if Murat.input.importLocation ~= 1
    
    eventStation                        =   zeros(lengthData,6);
    
else
    
    eventStation                        =   Murat.input.Locations;
    dist_xdeg_even                      =   eventStation(:,2)-origin(2);
    dist_ydeg_even                      =   eventStation(:,1)-origin(1);
    eventStation(:,3)                   =   -eventStation(:,3);
    
    dist_xdeg_station                   =   eventStation(:,5)-origin(2);
    dist_ydeg_station                   =   eventStation(:,4)-origin(1);
    
    %Transforms in meters
    eventStation(:,1)                   =   deg2km(dist_xdeg_even)*1000;
    eventStation(:,2)                   =   deg2km(dist_ydeg_even)*1000;
    eventStation(:,4)                   =   deg2km(dist_xdeg_station)*1000;
    eventStation(:,5)                   =   deg2km(dist_ydeg_station)*1000;

end

%Set up variables to save
location                                =   zeros(lengthData,6); 
theoreticalTime                         =   zeros(lengthData,1); 
peakDelay                               =   zeros(lengthData,1); 
totalLengthRay                          =   zeros(lengthData,1);
inverseQc                               =   zeros(lengthData,1); 
uncertaintyQc                           =   zeros(lengthData,1); 
travelTimePOrS                          =   zeros(lengthData,1);
energyRatioBodyCoda                     =   zeros(lengthData,1); 
energyRatioCodaNoise                    =   zeros(lengthData,1);
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
    
    %Display every 200 waveforms
    if isequal(mod(i,1000),0)
        
        disp(['Waveform number is ', num2str(i)])
        
    end
    
    %% OPERATIONS ON WAVEFORM
    listSac_i                           =   listSac{i};
    [pktime_i,sp_i,hsp_i,SAChdr_i,...
        srate_i,tempis]                 =...
        Murat_envelope(Murat,listSac_i);
    
    % In case it has been precalculated from external files.
    eventStation_i                      =   eventStation(i,:);
    
    if isequal(eventStation_i,zeros(1,6))
    
        location_i                      =   Murat_location(Murat,SAChdr_i);
        
    else
        
        location_i                      =   eventStation_i;
    
    end
    
    [theoreticalTime_i, originTime_i,...
        tCoda_i, cursorPick_i, cursorPeakDelay_i, cursorCodaStart_i,...
        cursorCodaEnd_i]                =...
        Murat_times(Murat,tempis,pktime_i,location_i,srate_i);
    travelTimePOrS_i                    =   pktime_i - originTime_i;
    
    %% OPERATIONS TO MEASURE AND MODEL PEAK DELAYS + INVERSION MATRIX Q
    peakDelay_i                         =...
        Murat_peakDelay(sp_i,cursorPick_i,cursorPeakDelay_i,srate_i);
    
    calculateRays                       =   recognizeComponents(i,compon);
    
    if calculateRays
        
        [Apd_i, AQ_i, totalLengthRay_i, raysPlot_i,...
            rayCrossing_i]              =...
            Murat_rays(Murat,location_i);
        
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
                Murat_codaMatrix(Murat,0,tCoda_i,location_i);
            
            inversionMatrixQc(i,:)      =   AQc_i;
    
    end
                
    %% OPERATIONS TO MEASURE  Q
    [energyRatioBodyCoda_i,...
        energyRatioCodaNoise_i]         = Murat_body(Murat,...
        srate_i,hsp_i,cursorPick_i,cursorCodaStart_i,cursorCodaEnd_i);
 
    %% SAVING
    location(i,:)                       =   location_i;
    theoreticalTime(i,1)                =   theoreticalTime_i;
    travelTimePOrS(i,1)                 =   travelTimePOrS_i;
    peakDelay(i,1)                      =   peakDelay_i;
    inverseQc(i,1)                      =   inverseQc_i; 
    uncertaintyQc(i,1)                  =   uncertaintyQc_i; 
    energyRatioBodyCoda(i,1)            =   energyRatioBodyCoda_i; 
    energyRatioCodaNoise(i,1)           =   energyRatioCodaNoise_i;

end

%% Setting up the final data vectors and matrices with checks on values
Murat.data.locations                    =   location;
Murat.data.theoreticalTime              =   theoreticalTime;
Murat.data.travelTimePOrS               =   travelTimePOrS;
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


function calculateValue         =...
    recognizeComponents(index,components)
%LOGICAL to decide if forward model is necessary depending in waveform
%number (index) and number of components.

calculateValue                  =   isequal(components,1) ||...
    (isequal(components,2) && ~isequal(mod(index,2),0)) ...
        || (isequal(components,1) && ~isequal(mod(index,2),0) &&...
        ~isequal(mod(index,3),0));

