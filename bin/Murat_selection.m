function Murat                      =   Murat_selection(Murat)

% Selects input and data to select
comp                                =   Murat.input.components;
tresholdnoise                       =   Murat.input.tresholdNoise;
modv                                =   Murat.input.modv;
nonlinear                           =   Murat.input.nonLinear;
PorS                                =   Murat.input.POrS;
listaSac                            =   Murat.input.listSac;
maPD                                =   Murat.input.maximumPeakDelay;
miPD                                =   Murat.input.minimumPeakDelay;

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

dataL                               =   size(peakd,1);
dataFreq                            =   size(peakd,2);
modvL                               =   size(modv,1);

fitrobust                           =   zeros(2,dataFreq);
lpdelta                             =   zeros(dataL,dataFreq);
retain_pd                           =   false(dataL,dataFreq);
retain_Qm                           =   false(dataL,dataFreq);
retain_Q                            =   false(dataL,dataFreq);
ray_crosses_pd                      =   false(modvL,dataFreq);
ray_crosses_Qc                      =   false(modvL,dataFreq);
ray_crosses_Q                       =   false(modvL,dataFreq);

switch nonlinear
    case 0
        fT                          =   Murat.input.fitTresholdLinear;
    case 1
        fT                          =   Murat.input.fitTresholdNonLinear;
end
    

%% Setting up the data vector in case of 2- and 3-components data

luntot                              =   luntot(1:comp:dataL);
time0                               =   tPS(1:comp:dataL);
evestaz                             =   evestaz(1:comp:dataL,:);
raysplot                            =   raysplot(:,:,1:comp:dataL);
Ac_i                                =   Ac_i(1:comp:dataL,:);
Apd_i                               =   Apd_i(1:comp:dataL,:);
A_i                                 =   A_i(1:comp:dataL,:);

if comp >  1
    [peakd,Qm,RZZ,rapsp,rapspcn]	=...
    Murat_components(components,peakd,Qm,RZZ,rapsp,rapspcn);
end


%% Warns about problematic data and saves their names
[problemPD,problemQc,problemRZZ,problemQ,yesPD]...
                                    =...
            Murat_dataWarning(listaSac,nonlinear,tresholdnoise,...
            maPD,miPD,fT,peakd,Qm,RZZ,rapspcn);


%% Operations to decide weight of each data for the solution
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
for i = 1:dataFreq
    
    % Peak Delays
    l10pd_i                         =   l10pd(:,i);
    yesPD_i                         =   yesPD(:,i);
    [pab,lpdelta_i,retain_pd_i,ray_crosses_pd_i]...
                                    =...
            Murat_retainPeakDealay(t_phase,l10pd_i,yesPD_i,Apd_i);
    
    % Qc
    Qm_i                            =   Qm(:,i);
    RZZ_i                           =   RZZ(:,i);
    [retain_Qm_i,ray_crosses_Qc_i]  =...
    Murat_retainQc(nonlinear,fT,Qm_i,RZZ_i,Ac_i);
    
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

