%% Seismic attributesare selected and components are considered
function Murat                      =   Murat_selection(Murat)
% SELECTS inputs and data
components                          =   Murat.input.components;
tresholdnoise                       =   Murat.input.tresholdNoise;
modv                                =   Murat.input.modv;
PorS                                =   Murat.input.POrS;
listaSac                            =   Murat.input.listSac;
maPD                                =   Murat.input.maximumPeakDelay;
miPD                                =   Murat.input.minimumPeakDelay;
fT                                  =   Murat.input.fitTresholdLinear;

Apd_i                               =...
    Murat.data.inversionMatrixPeakDelay;
A_i                                 =   Murat.data.inversionMatrixQ;
Ac_i                                =   Murat.data.inversionMatrixQc;
peakd                               =   Murat.data.peakDelay;
luntot                              =   Murat.data.totalLengthRay;
tPS                                 =   Murat.data.theoreticalTime;
evestaz                             =   Murat.data.locationsDeg;
Qm                                  =   Murat.data.inverseQc;
RZZ                                 =   Murat.data.uncertaintyQc;
rapsp                               =   Murat.data.energyRatioBodyCoda;
rapspcn                             =   Murat.data.energyRatioCodaNoise;
raysplot                            =   Murat.data.raysPlot;
tCoda                               =   Murat.data.tCoda;

dataL                               =   size(peakd,1);
dataFreq                            =   size(peakd,2);
modvL                               =   size(modv,1);

fitrobust                           =   zeros(2,dataFreq);
ray_crosses_pd                      =   false(modvL,dataFreq);
ray_crosses_Qc                      =   false(modvL,dataFreq);
ray_crosses_Q                       =   false(modvL,dataFreq);

%%
% Warns about problematic data and saves their names and locations
[problemPD,problemQc,problemRZZ,problemQ,~,compMissing,flagWarning]...
                                    =...
            Murat_dataWarning(listaSac,tresholdnoise,...
            maPD,miPD,fT,peakd,Qm,RZZ,rapspcn,components,0);

%%
% Selects data in case of multiple components
luntot                              =   luntot(1:components:dataL);
time0                               =   tPS(1:components:dataL);
tCoda                               =   tCoda(1:components:dataL,:);
evestaz                             =   evestaz(1:components:dataL,:);
raysplot                            =   raysplot(:,:,1:components:dataL);
Ac_i                                =   Ac_i(1:components:dataL,:);
Apd_i                               =   Apd_i(1:components:dataL,:);
A_i                                 =   A_i(1:components:dataL,:);

if components >  1
    [peakd,Qm,RZZ,rapsp,rapspcn]    =...
    Murat_components(components,peakd,Qm,RZZ,...
    rapsp,rapspcn,compMissing);
end
[~,~,~,~,yesPD,~,~]                 =...
            Murat_dataWarning(listaSac,tresholdnoise,...
            maPD,miPD,fT,peakd,Qm,RZZ,rapspcn,components,flagWarning);

%%
% Operations to decide weight of each data for the solution
% Using Vp/Vs to map max of S waves in the case of P picking
vpvs                                =   sqrt(3);
l10pd                               =   log10(peakd);

if PorS == 2
    time0                           =   time0*vpvs;
    t_phase                         =   log10(time0);
    
elseif PorS == 3
    t_phase                         =   log10(time0);

end

%%
% Remove outliers and inversion parameters with little/no sensitivity and
% store the remaining indexes for later
dataLMoreComp                       =   size(peakd,1);
lpdelta                             =   zeros(dataLMoreComp,dataFreq);
retain_pd                           =   false(dataLMoreComp,dataFreq);
retain_Qm                           =   false(dataLMoreComp,dataFreq);
retain_Q                            =   false(dataLMoreComp,dataFreq);
for i = 1:dataFreq
    
    % Peak Delays
    yesPD_i                         =   yesPD(:,i);
    l10pd_i                         =   l10pd(:,i);
    [pab,lpdelta_i,retain_pd_i,ray_crosses_pd_i]...
                                    =...
            Murat_retainPeakDelay(t_phase,l10pd_i,yesPD_i,Apd_i);
    
    % Qc
    Qm_i                            =   Qm(:,i);
    RZZ_i                           =   RZZ(:,i);
    [retain_Qm_i,ray_crosses_Qc_i]  =...
    Murat_retainQc(fT,Qm_i,RZZ_i,Ac_i);
    
    % Coda-normalization
    retain_Q_i                      =   rapspcn(:,i)>=tresholdnoise;
    A_retain_Q_i                    =   A_i(retain_Q_i,:);
    ray_crosses_Q_i                 =   sum(A_retain_Q_i)~=0;
    
    fitrobust(:,i)                  =   pab;
    lpdelta(:,i)                    =   lpdelta_i;
    retain_pd(:,i)                  =   retain_pd_i;
    ray_crosses_pd(:,i)             =   ray_crosses_pd_i;
    retain_Qm(:,i)                  =   retain_Qm_i;
    ray_crosses_Qc(:,i)             =   ray_crosses_Qc_i;
    retain_Q(:,i)                   =   retain_Q_i;
    ray_crosses_Q(:,i)              =   ray_crosses_Q_i;
   
end

Murat.data.peakDelay                =   peakd;
Murat.data.totalLengthRay           =   luntot;
Murat.data.tCoda                    =   tCoda;
Murat.data.locationsDeg             =   evestaz;
Murat.data.inverseQc                =   Qm;
Murat.data.uncertaintyQc            =   RZZ;
Murat.data.problemPD                =   problemPD;
Murat.data.problemQc                =   problemQc;
Murat.data.problemRZZ               =   problemRZZ;
Murat.data.problemQ                 =   problemQ;
Murat.data.energyRatioBodyCoda      =   rapsp;
Murat.data.energyRatioCodaNoise     =   rapspcn;
Murat.data.raysPlot                 =   raysplot;
Murat.data.variationPeakDelay       =   lpdelta;
Murat.data.travelTime               =   time0;
Murat.data.fitrobust                =   fitrobust;
Murat.data.retainPeakDelay          =   retain_pd;
Murat.data.retainQc                 =   retain_Qm;
Murat.data.retainQ                  =   retain_Q;
Murat.data.raysPeakDelay            =   ray_crosses_pd;
Murat.data.raysQc                   =   ray_crosses_Qc;
Murat.data.raysQ                    =   ray_crosses_Q;
Murat.data.inversionMatrixPeakDelay =   Apd_i;
Murat.data.inversionMatrixQc        =   Ac_i;
Murat.data.inversionMatrixQ         =   A_i;

end
