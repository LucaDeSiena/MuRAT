%% Peak-delay, Qc and Q TOMOGRAPHIC INVERSIONS
function Murat                      =   Murat_inversion(Murat)
%%
% Importing all the necessary inputs and data for plotting
FLabel                              =   Murat.input.label;
fformat                             =   Murat.input.format;
outputLCurve                        =   Murat.input.lCurve;
tWm                                 =   Murat.input.codaWindow;
cf                                  =   Murat.input.centralFrequency;
sped                                =   Murat.input.spectralDecay;
nxc                                 =   Murat.input.gridLong;
nyc                                 =   Murat.input.gridLat;
nzc                                 =   Murat.input.gridZ;
sizea                               =   Murat.input.sizeCheck;
latt                                =   Murat.input.lowCheck;
hatt                                =   Murat.input.highCheck;
modv                                =   Murat.input.modv;
spike_o                             =   Murat.input.spikeLocationOrigin;
spike_e                             =   Murat.input.spikeLocationEnd;
spike_v                             =   Murat.input.spikeValue;
x                                   =   Murat.input.x;
y                                   =   Murat.input.y;
z                                   =   Murat.input.z;
QcM                                 =   Murat.input.QcMeasurement;
inversionMethod                     =   Murat.input.inversionMethod;
lCurveQc                            =   Murat.input.lCurveQc;
lCurveQ                             =   Murat.input.lCurveQ;

Apd_i                               =...
    Murat.data.inversionMatrixPeakDelay;
Ac_i                                =   Murat.data.inversionMatrixQc;
A_i                                 =   Murat.data.inversionMatrixQ;
luntot                              =   Murat.data.totalLengthRay;
time0                               =   Murat.data.travelTime;
Qm                                  =   Murat.data.inverseQc;
RZZ                                 =   Murat.data.uncertaintyQc;
lpdelta                             =   Murat.data.variationPeakDelay;
rapsp                               =   Murat.data.energyRatioBodyCoda;
retain_pd                           =   Murat.data.retainPeakDelay;
retain_Qc                           =   Murat.data.retainQc;
retain_Q                            =   Murat.data.retainQ;
ray_crosses_pd                      =   Murat.data.raysPeakDelay;
ray_crosses_Qc                      =   Murat.data.raysQc;
ray_crosses_Q                       =   Murat.data.raysQ;
tCoda                               =   Murat.data.tCoda;
FPath                               =   './';

lMF                                 =   size(ray_crosses_pd);
modv_pd                             =   zeros(lMF(1),5,lMF(2));
modv_Qc                             =   zeros(lMF(1),10,lMF(2));
modv_Q                              =   zeros(lMF(1),10,lMF(2));
const_Qc                            =   zeros(size(rapsp));
residualQ                           =   zeros(1,lMF(2));
residualQc                          =   zeros(1,lMF(2));

%%
% Loops over all frequencies and parameter models
for k = 1:lMF(2)
    modv_pd(:,1:3,k)                =   modv;
    modv_Qc(:,1:3,k)                =   modv;
    modv_Q(:,1:3,k)                 =   modv;
    cf_k                            =   cf(k);
    fcName                          =   num2str(cf_k);
    if find(fcName == '.')
        fcName(fcName == '.')       =   '_';
    end
    
    %%
    % Peak delay regionalization
    rcpd_k                          =   ray_crosses_pd(:,k);
    rtpd_k                          =   retain_pd(:,k);
    Apd_k                           =	Apd_i(rtpd_k,rcpd_k);
    lpdelta_k                       =   lpdelta(rtpd_k,k);
    
    mpd                             =   sum(Apd_k.*lpdelta_k,1)';
    modv_pd(rcpd_k,4,k)             =   mpd;
    
    %%
    % Qc inversion
    rcQc_k                          =   ray_crosses_Qc(:,k);
    rtQc_k                          =   retain_Qc(:,k);
    Ac_k                            =   Ac_i(rtQc_k,rcQc_k);
    Qm_k                            =   Qm(rtQc_k,k);
    RZZ_k                           =   RZZ(rtQc_k,k);
    Wc                              =   Murat_weighting(RZZ_k,QcM);
    Gc                              =   Wc*Ac_k;
    FName                           =   ['L-curve_Qc_' fcName '_Hz'];
    
    bQm                             =   Wc*Qm_k;
    lCurveQc_k                      =   lCurveQc(k);
    
    if isequal(inversionMethod,'Tikhonov')
        [mtik0C,residualQc_k,LcQc,tik0_regC]...
                                    =...
           Murat_tikhonovQc(outputLCurve,Gc,bQm,lCurveQc_k);
        
        residualQc(1,k)             =   residualQc_k;
        modv_Qc(rcQc_k,4,k)         =   mtik0C;
        
    elseif isequal(inversionMethod,'Iterative')
        disp(['Qc L-curve and cost functions at ', num2str(cf_k), ' Hz.'])
        
        [LcQc, minimizeVectorQm,infoVectorQm,tik0_regC]...
                                    =...
           Murat_minimise(outputLCurve,Gc,bQm,lCurveQc_k,FName);
                                    
        residualQc(1,k)             =   min(infoVectorQm.Rnrm);
        modv_Qc(rcQc_k,4,k)         =   minimizeVectorQm;
        
    else
        error('Unknown inversion method.')
        
    end
    saveas(LcQc,fullfile(FPath, FLabel,'Rays_Checks',FName));
    close(LcQc)  
    
    %%
    % Q inversion
    rcQ_k                           =   ray_crosses_Q(:,k);
    rtQ_k                           =   retain_Q(:,k);
    A_k                             =   A_i(rtQ_k,rcQ_k);
    Q_k                             =   Qm(rtQ_k,k);
    luntot_k                        =   luntot(rtQ_k);
    time0_k                         =   time0(rtQ_k);
    rapsp_k                         =   rapsp(rtQ_k,k);    
    tCm                             =   tCoda(rtQ_k,k);
    
    [d1, const_Qc_k, ~, ~]          =   Murat_lsqlinQmean(tCm,tWm,Q_k,...
                                    cf_k,sped,luntot_k,time0_k,rapsp_k);
    const_Qc(rtQ_k,k)               =   const_Qc_k;
    
    lCurveQ_k                       =   lCurveQ(k);
    FName                           =   ['L-curve_Q_' fcName '_Hz'];
    if isequal(inversionMethod,'Tikhonov')
        [mtik0,residualQ_k,LcCN,tik0_reg]...
                                    =...
           Murat_tikhonovQ(outputLCurve,A_k,d1,lCurveQ_k,1);

        residualQ(:,k)              =   residualQ_k;
        modv_Q(rcQ_k,4,k)           =   mtik0;
        
    elseif isequal(inversionMethod,'Iterative')
        disp(['Q L-curve and cost functions at ', num2str(cf_k), ' Hz.'])
        [LcCN, minimizeVectorQ,infoVectorQ,tik0_reg]...
                                    =...
           Murat_minimise(outputLCurve,A_k,d1,lCurveQ_k,FName);
                                    
        residualQ(1,k)              =   min(infoVectorQ.Rnrm);
        modv_Q(rcQ_k,4,k)           =   minimizeVectorQ;
        
    end
    
    saveas(LcCN,fullfile(FPath, FLabel,'Rays_Checks',FName),fformat);
    close(LcCN)
    
    %% Checkerboards and spike inputs and checkerboard inversion
    % Qc
    siz                             =   [nxc nyc nzc];
    I                               =   checkerBoard3D(siz,sizea);
    [checkInput,spikeInput]         =...
                               Murat_inputTesting(I,spike_o,spike_e,x,y,z);
    
    modv_Qc(checkInput==1,6,k)      =   latt;
    modv_Qc(checkInput==0,6,k)      =   hatt;
    modv_Qc(:,8,k)                  =   mean(Qm_k);
    modv_Qc(spikeInput,8,k)         =   spike_v;
    Qc_ch                           =   modv_Qc(rcQc_k,6,k);
    re_checkQc                           =   Gc*Qc_ch;
    
    modv_Qc(rcQc_k,7,k)             =...
        Murat_outputTesting(Gc,re_checkQc,tik0_regC,inversionMethod);
    %%
    % Q
    modv_Q(:,6:8,k)                 =   modv_Qc(:,6:8,k);
    Q_ch                            =   modv_Q(rcQ_k,6,k);
    re_Q                            =   A_k*Q_ch;
    
    modv_Q(rcQ_k,7,k)               =...
        Murat_outputTesting(A_k,re_Q,tik0_reg,inversionMethod);
    
    %%
    % Inverting spike for Qc and Q at user discretion
    if ~isempty(spike_o)
        Qc_sp                       =   modv_Qc(rcQc_k,8,k);
        re_spikeQc                       =   Gc*Qc_sp;
        
        modv_Qc(rcQc_k,9,k)         =...
            Murat_outputTesting(Gc,re_spikeQc,tik0_regC,inversionMethod);
        
        Q_sp                        =   modv_Q(rcQ_k,8,k);
        re_spikeQ                        =   A_k*Q_sp;
        
        modv_Q(rcQ_k,9,k)           =...
            Murat_outputTesting(A_k,re_spikeQ,tik0_reg,inversionMethod);
        
    end
    
    %%
    % Save peak-delay, Qc, Q
    modv_pd_k                       =   modv_pd(:,:,k);
    FName                           =   ['peakdelay_' fcName '_Hz.txt'];
    save(fullfile(FPath, FLabel, 'TXT', FName), 'modv_pd_k','-ascii');
    
    modv_Qc_k                       =   modv_Qc(:,:,k);
    FName                           =   ['Qc_' fcName '_Hz.txt'];
    save(fullfile(FPath, FLabel, 'TXT', FName), 'modv_Qc_k','-ascii');
    
    modv_Q_k                        =   modv_Q(:,:,k);
    FName                           =   ['Q_' fcName '_Hz.txt'];
    save(fullfile(FPath, FLabel, 'TXT', FName), 'modv_Q_k','-ascii');
    
end
%%
% Save in Murat
Murat.data.residualQc               =   residualQc;
Murat.data.const_Qc                 =   const_Qc;
Murat.data.residualQ                =   residualQ;
Murat.data.modvPeakDelay            =   modv_pd;
Murat.data.modvQc                   =   modv_Qc;
Murat.data.modvQ                    =   modv_Q;
