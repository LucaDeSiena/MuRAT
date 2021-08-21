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

if outputLCurve == 0
    lCurveQc                        =   Murat.input.lCurveQc;
    lCurveQ                         =   Murat.input.lCurveQ;
else
    lCurveQc                        =   [];
    lCurveQ                         =   [];
end

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
    modv_pd(:,1:4,k)                =   modv;
    modv_Qc(:,1:4,k)                =   modv;
    modv_Q(:,1:4,k)                 =   modv;
    
    %%
    % Peak delay regionalization
    rcpd                            =   ray_crosses_pd(:,k);
    rtpd                            =   retain_pd(:,k);
    Apd                             =	Apd_i(rtpd,rcpd);
    lApd                            =   size(Apd,2);
    mpd                             =   zeros(lApd,1);
    cfk                             =   cf(k);
    fcName                          =   num2str(cfk);
    
    for j = 1:lApd
        
        mpd(j,1)                    =...
            sum(Apd(:,j).*lpdelta(rtpd,k))/sum(Apd(:,j));
        
    end
    
    modv_pd(rcpd,5,k)               =   mpd;
    
    %%
    % Qc inversion
    rcQc                            =   ray_crosses_Qc(:,k);
    rtQc                            =   retain_Qc(:,k);
    Ac                              =   Ac_i(rtQc,rcQc);
    Qm_k                            =   Qm(rtQc,k);
    RZZ_k                           =   RZZ(rtQc,k);
    Wc                              =   Murat_weighting(RZZ_k,QcM);
    Gc                              =   Wc*Ac;
    [Uc,Sc,Vc]                      =   svd(Gc);
    
    [mtik0C,residualQc_k,LcQc,tik0_regC]...
                                    =...
        Murat_tikhonovQc(outputLCurve,Qm_k,Wc,Gc,lCurveQc);
    
    FName                           =   ['L-curve_Qc_' fcName '_Hz'];
    saveas(LcQc,fullfile(FPath, FLabel,'Rays_Checks',FName),fformat);
    close(LcQc)

    residualQc(1,k)                 =   residualQc_k;
    modv_Qc(rcQc,5,k)               =   mtik0C;
    
    %%
    % Q inversion
    rcQ                             =   ray_crosses_Q(:,k);
    rtQ                             =   retain_Q(:,k);
    A                               =   A_i(rtQ,rcQ);
    Q_k                             =   Qm(rtQ,k);
    
    tCm                             =   mean(tCoda(:,k));
    const_Qc_k                      =...
        (tCm+tWm/2)^-sped.*exp(-Q_k*2*pi*cfk*(tCm+tWm/2));
    rapsp_k                         =   rapsp(rtQ,k);    
    [U,S,V]                         =   svd(A);
    
    [mtik0,residualQ_k,LcCN,tik0_reg,~,~]...
                                    =   Murat_tikhonovQ(cfk,rtQ,...
        outputLCurve,rapsp_k,const_Qc_k,luntot,time0,A,lCurveQ,1);

    FName                           =   ['L-curve_Q_' fcName '_Hz'];
    saveas(LcCN,fullfile(FPath, FLabel,'Rays_Checks',FName),fformat);
    close(LcCN)
    modv_Q(rcQ,5,k)                 =   mtik0;
    const_Qc(rtQ,k)                 =   const_Qc_k;
    residualQ(:,k)                  =   residualQ_k;
    
    %%
    % Checkerboards and spike inputs and checkerboard inversion
    siz                             =   [nxc nyc nzc];
    I                               =   checkerBoard3D(siz,sizea);
    [checkInput,spikeInput]         =...
    Murat_inputTesting(I,spike_o,spike_e,x,y,z);
    
    modv_Qc(checkInput==1,6,k)      =   latt;
    modv_Qc(checkInput==0,6,k)      =   hatt;
    modv_Qc(:,8,k)                  =   mean(Qm_k);
    modv_Qc(spikeInput,8,k)         =   spike_v;
    modv_Q(:,6:8,k)                 =   modv_Qc(:,6:8,k);
    Qc_ch                           =   modv_Qc(rcQc,6,k);
    re_Qc                           =   Gc*Qc_ch;
    mcheck_c                        =...
        tikhonov(Uc,diag(Sc),Vc,re_Qc,tik0_regC);
    modv_Qc(rcQc,7,k)               =   mcheck_c;
    Q_ch                            =   modv_Q(rcQ,6,k);
    re_Q                            =   A*Q_ch;
    mcheck                          =...
        tikhonov(U,diag(S),V,re_Q,tik0_reg);
    modv_Q(rcQ,7,k)                 =   mcheck;
    
    %%
    % Inverting spike for Qc and Q at user discretion
    if ~isempty(spike_o)
        Qc_sp                       =   modv_Qc(rcQc,8,k);
        re_Qc                       =   Gc*Qc_sp;
        mspike_c                    =   tikhonov(Uc,diag(Sc),Vc,re_Qc,tik0_regC);
        modv_Qc(rcQc,9,k)           =   mspike_c;
        
        Q_sp                        =   modv_Q(rcQ,8,k);
        re_Q                        =   A*Q_sp;
        mspike                      =   tikhonov(U,diag(S),V,re_Q,tik0_reg);
        modv_Q(rcQ,9,k)             =   mspike;
    end
    
    %%
    % Diagonal of resolution matrix
    sSc                             =   size(Sc);
    fil_reg                         =   fil_fac(diag(Sc),tik0_regC);
    if sSc(2) > sSc(1)
        fil_reg(end+1:sSc(2),1)     =   0;
    end
    Rc                              =   Vc*diag(fil_reg)*Vc';
    dRc                             =   diag(Rc);
    modv_Qc(rcQc,10,k)              =   dRc;
    
    sS                              =   size(S);
    fil_reg                         =   fil_fac(diag(S),tik0_reg);
    if sS(2) > sS(1)
        fil_reg(end+1:sS(2),1)      =   0;
    end
    R                               =   V*diag(fil_reg)*V';
    dR                              =   diag(R);
    modv_Q(rcQ,10,k)                =   dR;
    
    %%
    % Save peak-delay, Qc, Q
    modv_pd_k                       =   modv_pd(:,:,k);
    FName                           =   ['peakdelay_' fcName '_Hz.txt'];
    save(fullfile(FPath, FLabel, 'TXT', FName), 'modv_pd_k','-ascii');
    
    modv_Qc_k                       =   modv_pd(:,:,k);
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
