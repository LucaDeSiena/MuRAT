function Murat                      = Murat_selection(Murat)

% Selects input and data to select
compon                              =   Murat.input.components;
tresholdnoise                       =   Murat.input.tresholdNoise;
modv                                =   Murat.input.modv;
nonlinear                           =   Murat.input.nonLinear;
PorS                                =   Murat.input.POrS;

Apd                                 =...
    Murat.data.inversionMatrixPeakDelay;
A                                   =   Murat.data.inversionMatrixQ;
Ac                                  =   Murat.data.inversionMatrixQc;
peakd                               =   Murat.data.peakDelay;
luntot                              =   Murat.data.totalLengthRay;
tPS                                 =   Murat.data.theoreticalTime;
evestaz                             =   Murat.data.locationsDeg;
Qm                                  =   Murat.data.inverseQc;
RZZ                                 =   Murat.data.uncertaintyQc;
rapsp                               =   Murat.data.energyRatioBodyCoda;
rapspcn                             =   Murat.data.energyRatioCodaNoise;
raysplot                            =   Murat.data.raysPlot;

%Find the amount of Qc unavailable in the data
noQm                                =   sum(Qm==0)/length(Qm)*100;

%Conditions for the linear and non-linear inversions are opposite
if nonlinear == 0
    RZZ2                            =   Murat.input.fitTresholdLinear;
    noRZZ                           =   sum(RZZ<=RZZ2)/length(RZZ)*100;
elseif nonlinear == 1
    RZZ2                            =   Murat.input.fitTresholdNonLinear;
    noRZZ                           =   sum(RZZ>=RZZ2)/length(RZZ)*100;
end

if noQm > 20
    Qless                           =   ...
        [num2str(noQm),'% of your Q are =0'];
    warning(Qless)
end

if noRZZ>20
    Rless                           =   [num2str(noRZZ),...
        '% of your coerrelation coefficients are lower than treshold'];
    warning(Rless)
end

%% Setting up the data vector in case of 2- and 3-components data

%Only one waveform component
dataL                               =   length(peakd);
luntot                              =   luntot(1:compon:dataL);
time0                               =   tPS(1:compon:dataL);
evestaz                             =   evestaz(1:compon:dataL,:);
raysplot                            =   raysplot(:,:,1:compon:dataL);
Ac                                  =   Ac(1:compon:dataL,:);
Apd                                 =   Apd(1:compon:dataL,:);
A                                   =   A(1:compon:dataL,:);

% 2 components (typically [WE, SN]) - average the data
if compon ==  2
    
    peakd                           =...
        (peakd(1:comp:dataL-1) + peakd(2:comp:dataL))/2;
    Qm                              =...
        (Qm(1:comp:dataL-1) + Qm(2:comp:dataL))/2;
    RZZ                             =...
        (RZZ(1:comp:dataL-1) + RZZ(2:comp:dataL))/2;
    rapsp                           =...
        (rapsp(1:comp:dataL-1) + rapsp(2:comp:dataL))/2;
    rapspcn                         =...
        (rapspcn(1:comp:dataL-1) + rapspcn(2:comp:dataL))/2;
    
% 3 components (WE, SN, Z) - average the H data first
elseif compon == 3
    
    %Averaging the horizontals before the vertical
    peakd                           =...
        ((peakd(1:comp:dataL-2) + peakd(2:comp:dataL-1))/2 +....
        peakd(3:comp:dataL))/2;
    Qm                              =...
        ((Qm(1:comp:dataL-2) + Qm(2:comp:dataL-1))/2 +....
        Qm(3:comp:dataL))/2;
    RZZ                             =...
        ((RZZ(1:comp:dataL-2) + RZZ(2:comp:dataL-1))/2 +....
        RZZ(3:comp:dataL))/2;
    rapsp                           =...
        ((rapsp(1:comp:dataL-2) + rapsp(2:comp:dataL-1))/2 +....
        rapsp(3:comp:dataL))/2;
    rapspcn                         =...
        ((rapspcn(1:comp:dataL-2) + rapspcn(2:comp:dataL-1))/2 +....
        rapspcn(3:comp:dataL))/2;
    
end

%% Peak-delay
% setting up the trehold for acceptance - minimum peak delay 
yes_pd                              =...
    peakd > Murat.input.minimumPeakDelay;

%Using Vp/Vs to map max of S waves in case of P
vpvs                                =   sqrt(3);
l10pd                               =   log10(peakd);
%Case for P waves
if PorS == 2
    time0                           =   time0*vpvs;
    t_phase                         =   log10(time0);
    
    % Case for S waves
elseif PorS == 3
    t_phase                         =   log10(time0);
end

% Fitting peak delay values
fitrobust                           =...
    fit(t_phase,l10pd,'poly1','Robust','on');

% Coefficients of the linear relationship
pab                                 =   [fitrobust.p1 fitrobust.p2];

% Obtaining theoretical values and residuals
l10pdt                              =   polyval(pab,t_phase);
lpdelta                             =   l10pd-l10pdt;

% To remove outliers of peak delay
I                                   =   abs(lpdelta) >= 2*std(lpdelta);
outlierspd                          =...
    excludedata(t_phase,lpdelta,'indices',I);
retain_pd                           =   yes_pd & ~outlierspd;

%Remove anomalous Qc remembering that there are two different types of
%inversion

if nonlinear == 0
    retain_Qm                       =   find(Qm>0 & RZZ>RZZ2);
    mQm                             =   mean(Qm(retain_Qm));
    outliersc                       =   Qm > mQm+2*std(Qm(retain_Qm));
    retain_Qm                       =...
        find(Qm>0 & RZZ>RZZ2 & outliersc==0);
    
elseif nonlinear == 1
    retain_Qm                       =   find(Qm>0 & RZZ<RZZ2);
    mQm                             =   mean(Qm(retain_Qm));
    outliersc                       =   Qm > mQm+2*std(Qm(retain_Qm));
    retain_Qm                       =...
        find(Qm>0 & RZZ<RZZ2 & outliersc==0);
end

%% Remove outliers and inversion parameters with little/no sensitivity
% Peak delay data
Apd                                 =   Apd(retain_pd,:);%no outliers

% Peak delay inversion matrix and relative coords
s_pd                                =   sum(Apd); %find zero sensitivity
ray_crosses_pd                      =   s_pd~=0;
Apd                                 =   Apd(:,ray_crosses_pd);
modv_pd                             =   modv(:,1:3);
modv_pd(ray_crosses_pd,5)           =   1;

% Q data
retain_Q                            =   rapspcn>=tresholdnoise; %min coda-noise
A                                   =   A(retain_Q,:);

% Q inversion matrix and relative coords
s_Q                                 =   sum(A); %find zero sensitivity
ray_crosses_Q                       =   s_Q~=0;
A                                   =   A(:,ray_crosses_Q);
modv_Q                              =   modv(:,1:3);
modv_Q (ray_crosses_Q,5)            =   1;

% Qc data
Ac                                  =   Ac(retain_Qm,:);

% Qc inversion matrix
Ac(Ac<10^(-3))                      =   0;
s_Qc                                =   sum(Ac); %same for A
K_crosses                           =   s_Qc~=0;
Ac                                  =   Ac(:,K_crosses);
modv_Qc                             =   modv(:,1:3);
modv_Qc(K_crosses,5)                =   1;

Murat.data.inversionMatrixPeakDelay =   Apd;
Murat.data.inversionMatrixQ         =   A;
Murat.data.inversionMatrixQc        =   Ac;
Murat.data.peakDelay                =   peakd;
Murat.data.totalLengthRay           =   luntot;
Murat.data.locationsDeg             =   evestaz;
Murat.data.inverseQc                =   Qm;
Murat.data.uncertaintyQc            =   RZZ;
Murat.data.energyRatioBodyCoda      =   rapsp;
Murat.data.energyRatioCodaNoise     =   rapspcn;
Murat.data.raysPlot                 =   raysplot;
Murat.data.plotRays                 =   raysplot;

Murat.data.variationPeakDelay       =   lpdelta;
Murat.data.travelTime               =   time0;
Murat.data.fitrobust                =   fitrobust;
Murat.data.retainPeakDelay          =   retain_pd;
Murat.data.retainQc                 =   retain_Qm;
Murat.data.retainQ                  =   retain_Q;
Murat.data.modvPeakDelay            =   modv_pd;
Murat.data.modvQ                    =   modv_Q;
Murat.data.modvQc                   =   modv_Qc;

