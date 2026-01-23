%% Seismic attributesare selected and components are considered
function Murat  =   Murat_selection(Murat)
% SELECTS inputs and data
s               =   Murat.input;

FLabel          =   s.label;
components      =   s.components;
tresholdnoise   =   s.tresholdNoise;
modv            =   s.modv;
PorS            =   s.POrS;
listaSac        =   s.listSac;
maPD            =   s.maximumPeakDelay;
miPD            =   s.minimumPeakDelay;
fT              =   s.fitTresholdLinear;
stDevPD         =   s.stDevPD;
origin          =   s.origin;
ending          =   s.end;
x               =   s.x;
y               =   s.y;
z               =   s.z;
QcM             =   s.QcMeasurement;
cf              =   s.centralFrequency;

Declust         =   s.declustering;

Apd_i           =   Murat.PD.inversionMatrixPeakDelay;
A_i             =   Murat.Q.inversionMatrixQ;
Ac_i            =   Murat.Qc.inversionMatrixQc;
peakd           =   Murat.PD.peakDelay;
luntot          =   Murat.rays.totalLengthRay;
tPS             =   Murat.rays.theoreticalTime;
locationsDeg    =   Murat.rays.locationsDeg;
Qm              =   Murat.Qc.inverseQc;
RZZ             =   Murat.Qc.uncertaintyQc;
rapsp           =   Murat.Q.energyRatioBodyCoda;
rapspcn         =   Murat.Q.energyRatioCodaNoise;
raysplot        =   Murat.rays.raysPlot;
tCoda           =   Murat.Qc.tCoda;

dataL           =   size(peakd,1);
dataFreq        =   size(peakd,2);
modvL           =   size(modv,1);

fitrobust       =   zeros(2,dataFreq);
ray_crosses_pd  =   false(modvL,dataFreq);
ray_crosses_Qc  =   false(modvL,dataFreq);
ray_crosses_Q   =   false(modvL,dataFreq);

if ~isempty(Declust)
    disp(['Waveforms used before Qc and Q declustering: ',...
            num2str(numel(listaSac))])
    s.locDegOriginal=   locationsDeg;
    % prepare folder path builder for clustering outputs (used only if Declust)
    storeFolder = 'Tests/Clustering';
    baseFolder  = fullfile('./', FLabel, storeFolder);

    
end
        
%%
% Warns about problematic data and saves their names and locations
[problemPD,problemQc,problemRZZ,problemQ,~,compMissing,flagWarning]=...
            Murat_dataWarning(listaSac,tresholdnoise,...
            maPD,miPD,fT,peakd,Qm,RZZ,rapspcn,rapsp,components,0);

%%
% Selects data in case of multiple components
idxComp         =   1:components:dataL;

luntot          =   luntot(idxComp);
time0           =   tPS(idxComp);
tCoda           =   tCoda(idxComp,:);
locationsDeg    =   locationsDeg(idxComp,:);
raysplot        =   raysplot(:,:,idxComp);
Ac_i            =   Ac_i(idxComp,:,:);
Apd_i           =   Apd_i(idxComp,:);
A_i             =   A_i(idxComp,:);

if components >  1
    [peakd,Qm,RZZ,rapsp,rapspcn] =   Murat_components(components,peakd,...
        Qm,RZZ,rapsp,rapspcn,compMissing);
end
[~,~,~,~,yesPD,~,~] =    Murat_dataWarning(listaSac,tresholdnoise,maPD,...
    miPD,fT,peakd,Qm,RZZ,rapspcn,rapsp,components,flagWarning);

%%
% Operations to decide weight of each data for the solution
% Using Vp/Vs to map max of S waves in the case of P picking
vpvs            =   sqrt(3);
l10pd           =   log10(peakd);

if PorS == 2
    time0       =   time0*vpvs;
    t_phase     =   log10(time0);
    
elseif PorS == 3
    t_phase     =   log10(time0);

end

%%
% Remove outliers and inversion parameters with little/no sensitivity and
% store the remaining indexes for later
dataLMoreComp   =   size(peakd,1);
lpdelta         =   zeros(dataLMoreComp,dataFreq);
retain_pd       =   false(dataLMoreComp,dataFreq);
retain_Qm       =   false(dataLMoreComp,dataFreq);
retain_Q        =   false(dataLMoreComp,dataFreq);

for i = 1:dataFreq
    cf_k                            =   cf(i);
    fcName                          =   num2str(cf_k);
    if find(fcName == '.')
        fcName(fcName == '.')       =   '_';
    end
    
    % Peak Delays
    yesPD_i     =   yesPD(:,i);
    l10pd_i     =   l10pd(:,i);
    [pab,lpdelta_i,retain_pd_i,ray_crosses_pd_i]=...
            Murat_retainPeakDelay(t_phase,l10pd_i,yesPD_i,Apd_i,stDevPD);
    
    % Qc
    Qm_i        =   Qm(:,i);
    RZZ_i       =   RZZ(:,i);
    retain_Qc_i =   Murat_retainQc(fT,Qm_i,RZZ_i,QcM);

    if ~isempty(Declust)
        
        indKeepQc       =   Murat_declustering(origin,ending,x,y,z,...
            QcM,RZZ_i,locationsDeg,Declust);
        retain_Qc_i     =   retain_Qc_i & indKeepQc;

        FName_Cluster   =   ['ClusteringQc_' fcName '_Hz'];
        
        clustering      =   Murat_imageDeclustering(...
        locationsDeg,locationsDeg(indKeepQc,:),origin,ending,FName_Cluster);
        pathFolder = fullfile(baseFolder, FName_Cluster);
        saveFigureAsImage(pathFolder);
        close(clustering)

    end
    Ac_retain_Qc_i      =   Ac_i(retain_Qc_i,:,i);
    ray_crosses_Qc_i    =   sum(Ac_retain_Qc_i) > 1e-19;
    
    % Coda-normalization
    retainQCodaNoise    =   rapspcn(:,i) >= tresholdnoise;
    retainQNaN          =   ~isnan(Qm_i);
    retainQPSCoda       =   rapsp(:,i) >= tresholdnoise;
    retain_Q_i          =...
        retainQNaN & retainQCodaNoise & retainQPSCoda;
    
    if ~isempty(Declust)
        
        indKeepQ        =   Murat_declustering(origin,ending,x,y,z,...
            QcM,1./rapspcn(:,i),locationsDeg,Declust);
        retain_Q_i      =   retain_Q_i & indKeepQ;
        
        FName_Cluster   =   ['ClusteringQ_' fcName '_Hz'];
        clustering      =   Murat_imageDeclustering(...
        locationsDeg,locationsDeg(indKeepQ,:),origin,ending,FName_Cluster);
        pathFolder = fullfile(baseFolder, FName_Cluster);
        saveFigureAsImage(pathFolder);
        close(clustering)

    end
    
    A_retain_Q_i        =   A_i(retain_Q_i,:);
    ray_crosses_Q_i     =   sum(A_retain_Q_i) ~= 0;
    
    fitrobust(:,i)      =   pab;
    lpdelta(:,i)        =   lpdelta_i;
    retain_pd(:,i)      =   retain_pd_i;
    ray_crosses_pd(:,i) =   ray_crosses_pd_i;
    retain_Qm(:,i)      =   retain_Qc_i;
    ray_crosses_Qc(:,i) =   ray_crosses_Qc_i;
    retain_Q(:,i)       =   retain_Q_i;
    ray_crosses_Q(:,i)  =   ray_crosses_Q_i;
   
end

disp(['Waveforms used after Qc declustering: ', num2str(sum(retain_Qm))])
disp(['Waveforms used after Q declustering: ', num2str(sum(retain_Q))])

Murat.rays.locationsDeg     =   locationsDeg;
Murat.rays.totalLengthRay   =   luntot;
Murat.rays.raysPlot         =   raysplot;
Murat.rays.travelTime       =   time0;

Murat.PD.peakDelay          =   peakd;
Murat.PD.problemPD          =   problemPD;
Murat.PD.variationPeakDelay =   lpdelta;
Murat.PD.fitrobust          =   fitrobust;
Murat.PD.retainPeakDelay    =   retain_pd;
Murat.PD.raysPeakDelay      =   ray_crosses_pd;
Murat.PD.inversionMatrixPeakDelay   =   Apd_i;

Murat.Qc.tCoda              =   tCoda;
Murat.Qc.inverseQc          =   Qm;
Murat.Qc.uncertaintyQc      =   RZZ;
Murat.Qc.problemQc          =   problemQc;
Murat.Qc.problemRZZ         =   problemRZZ;
Murat.Qc.raysQc             =   ray_crosses_Qc;
Murat.Qc.retainQc           =   retain_Qm;
Murat.Qc.inversionMatrixQc  =   Ac_i;

Murat.Q.problemQ            =   problemQ;
Murat.Q.energyRatioBodyCoda =   rapsp;
Murat.Q.energyRatioCodaNoise=   rapspcn;
Murat.Q.retainQ             =   retain_Q;
Murat.Q.raysQ               =   ray_crosses_Q;
Murat.Q.inversionMatrixQ    =   A_i;

end
