function Murat                  =   Murat_inversion(Murat)
%  2D peak-delay, Qc and Q TOMOGRAPHIC INVERSIONS
% Using Tikhonov inversions and L curves

%INPUTS
FPath                           =   Murat.input.workingDirectory;
FLabel                          =   Murat.input.label;
fformat                         =   Murat.input.format;
outputLCurve                    =   Murat.input.lCurve;
tCm                             =   Murat.input.startLapseTime;
tWm                             =   Murat.input.codaWindow;
cf                              =   Murat.input.centralFrequency;
sped                            =   Murat.input.spectralDecay;
nxc                             =   Murat.input.gridLong;
nyc                             =   Murat.input.gridLat;
nzc                             =   Murat.input.gridZ;
sizea                           =   Murat.input.sizeCheck;
latt                            =   Murat.input.lowCheck;
hatt                            =   Murat.input.highCheck;
nonlinear                       =   Murat.input.nonLinear;
spike_o                         =   Murat.input.spikeLocationOrigin;
spike_e                         =   Murat.input.spikeLocationEnd;
spike_v                         =   Murat.input.spikeValue;
x                               =   Murat.input.x;
y                               =   Murat.input.y;
z                               =   Murat.input.z;

Apd                             =   Murat.data.inversionMatrixPeakDelay;
Ac                              =   Murat.data.inversionMatrixQc;
A                               =   Murat.data.inversionMatrixQ;
luntot                          =   Murat.data.totalLengthRay;
time0                           =   Murat.data.travelTime;
Qm                              =   Murat.data.inverseQc;
RZZ                             =   Murat.data.uncertaintyQc;
lpdelta                         =   Murat.data.variationPeakDelay;
rapsp                           =   Murat.data.energyRatioBodyCoda;
retain_pd                       =   Murat.data.retainPeakDelay;
retain_Qc                       =   Murat.data.retainQc;
retain_Q                        =   Murat.data.retainQ;
modv_pd                         =   Murat.data.modvPeakDelay;
modv_Qc                         =   Murat.data.modvQc;
modv_Q                          =   Murat.data.modvQ;


%% Peak delay mapping - weighted average
lApd                            =   size(Apd);
mpd                             =   zeros(lApd(2),1); %4th column peak delay

% For loop to sum and normalize
for j = 1:lApd(2)
    mpd(j,1)                    =...
        sum(Apd(:,j).*lpdelta(retain_pd))/sum(Apd(:,j));
end
%Assign to 4th column of what you save
modv_pd(modv_pd(:,5)==1,4)      =   mpd;

%% Qc inversion
% Weighting depends on RZZ
if nonlinear == 0
    
    W1                          =   RZZ(retain_Qc)<0.3;
    W2                          =   RZZ(retain_Qc)<0.5;
    W3                          =   RZZ(retain_Qc)<0.7;
    W4                          =   W3+W2+W1;
    Wc                          =   diag(1./(W4+1));% weights
    
    % For the nonlinear case there is no weighting
elseif nonlinear == 1
    
    Wc                          =   1;
    
end

%Weighted tikhonov
dcW                             =   Wc*Qm(retain_Qc);
Gc                              =   Wc*Ac;
[Uc,Sc,Vc]                      =   svd(Gc);

if outputLCurve == 1
    
    % Damped parameter is user defined
    LcQc                        =...
        figure('Name','L-curve Qc','NumberTitle','off');
    [rho,eta,reg_param]         =...
        l_curve_tikh_svd(Uc,diag(Sc),dcW,100);
    plot_lc(rho,eta,'-',1,reg_param)
%     l_curve(Uc,diag(Sc),dcW,'Tikh')
    % define from L-curve
    tik0_regC                   =...
        input('Your personal smoothing parameter for coda ');
    
    FName                       =   'Lc_Qc';
    saveas(LcQc,fullfile(FPath, FLabel, FName), fformat);
    close(LcQc)
    
    % picard plot
    PpQc                        =   figure('Name','Picard-plot Qc',...
        'NumberTitle','off','visible','off');
    picard(Uc,diag(Sc),dcW);
    FName                       =   'Picard_Qc';
    saveas(PpQc,fullfile(FPath, FLabel, FName), fformat);
else
    tik0_regC                   =   Murat.input.lCurveQc;
    
end
% invert
mtik0C                          =   tikhonov(Uc,diag(Sc),Vc,dcW,tik0_regC);
modv_Qc(modv_Qc(:,5)==1,4)      =   mtik0C;

%% Q inversion - Direct wave attenuation

% Set up the constant from start - the Qc along paths are needed
const_Qc                        =...
    (tCm+tWm/2)^sped.*exp(Qm(retain_Q)*2*pi*cf*(tCm+tWm/2));

% Data of the inverse problem
d0                              =   log(rapsp(retain_Q)./const_Qc)/2/pi/cf;

% Inversion matrix - geom. spre. + average Q
G                               =   -log(luntot(retain_Q))/pi/cf;
G(:,2)                          =   -time0(retain_Q);

% The 3 parameters are the constant, the geometrical spreading, and the
% average Q and they are contained in constQmean.
constQmean                      =   lsqlin(G,d0(:,1));% least square inversion
cova                            =...
    (G'*G)^(-1)*G'*cov(d0)*G*(G'*G)^(-1); %covariance matrix
er                              =   sqrt(diag(cova)); %error from the covariance matrix
constQmean(:,2)                 =   er;
%Data creation, removing the pre-calculated parameters
d1                              =...
    d0  + constQmean(1,1)*log(luntot(retain_Q))/pi/cf...
    + time0(retain_Q)*constQmean(2,1);


% NEW DATA VECTOR
% tikhonov inversion - by using the programs in HANSEN et al. 1994
[U,S,V]=svd(A);

if outputLCurve == 1
    
    %sets the smoothing parameter - always user defined
    LcCN                        =   figure('Name','L-curve coda-normalization','NumberTitle','off');
    l_curve(U,diag(S),d1,'Tikh')
    tik0_reg                    =   input('Your personal smoothing parameter ');
    FName                       =   'Lc_CN';
    saveas(LcCN,fullfile(FPath, FLabel, FName), fformat);
    close(LcCN)
    
    % picard plot
    PpCN                        =   figure('Name','Picard coda-norm.','NumberTitle','off',...
        'visible','off');
    picard(U,diag(S),d1);
    FName                       =   'Picard_CN';
    saveas(PpCN,fullfile(FPath, FLabel, FName), fformat);
    
else
    tik0_reg                    =   Murat.input.lCurveQ;
    
end
%results
mtik0=tikhonov(U,diag(S),V,d1,tik0_reg);
modv_Q(modv_Q(:,5)==1,4)        =   mtik0;

%% Testing - 3D checkerboard and spike inputs

% Start with checkerboard
siz                             =   [nxc nyc nzc];

% Function to create 3D pattern
I                               =   checkerBoard3D(siz,sizea);

% Unfold in checkerboard velocity model
index                           = 0;
for i=1:nxc
    for j=1:nyc
        for k=1:nzc
            index               = index+1;
            modv_Qc(index,6)    = I(i,j,k);
        end
    end
end

% Checkeroard inside Qc - input
modv_Qc(modv_Qc(:,6)==1,6)      =   latt;
modv_Qc(modv_Qc(:,6)==0,6)      =   hatt;


% Spike test
modv_Qc(:,8)                    =   mean(Qm);

if ~isempty(spike_o)
    r                           =   Murat_unfold(x',y',z');
    cond_in                     =...
        r(:,1)>spike_o(2) & r(:,1)<spike_e(2) &...
        r(:,2)>spike_o(1) & r(:,2)<spike_e(1) &...
        r(:,3)<spike_o(3) & r(:,3)>spike_e(3);
    
    modv_Qc(cond_in,8)          =   spike_v;
end

% Checkerboard/spike for Q - input is the same
modv_Q(:,6:8)                   =   modv_Qc(:,6:8);

%% Testing - 3D checkerboard and spike inversions
% invert checkerboard for Qc
Qc_ch                           =   modv_Qc(modv_Qc(:,5)==1,6);
re_Qc                           =   Gc*Qc_ch;
mcheck_c                        =   tikhonov(Uc,diag(Sc),Vc,re_Qc,tik0_regC);
modv_Qc(modv_Qc(:,5)==1,7)      =   mcheck_c; %output checkerboard

% invert checkerboard for Q
Q_ch                            =   modv_Q(modv_Q(:,5)==1,6);
re_Q                            =   A*Q_ch;
mcheck                          =   tikhonov(U,diag(S),V,re_Q,tik0_reg);
modv_Q(modv_Q(:,5)==1,7)        =   mcheck; %output checkerboard

if ~isempty(spike_o)
    % invert spike for Qc
    Qc_sp                       =   modv_Qc(modv_Qc(:,5)==1,8);
    re_Qc                       =   Gc*Qc_sp;
    mspike_c                    =   tikhonov(Uc,diag(Sc),Vc,re_Qc,tik0_regC);
    modv_Qc(modv_Qc(:,5)==1,9)  =   mspike_c; %output spike
    
    % invert spike for Qc
    Q_sp                        =   modv_Q(modv_Q(:,5)==1,8);
    re_Q                        =   A*Q_sp;
    mspike                      =   tikhonov(U,diag(S),V,re_Q,tik0_reg);
    modv_Q(modv_Q(:,5)==1,9)    =   mspike; %output spike
end

% WE ALSO OUTPUT THE DIAGONAL OF THE RESOLUTION MATRIX
% using the filter functions for Qc
sSc                             =   size(Sc);
fil_reg                         =   fil_fac(diag(Sc),tik0_regC);
if sSc(2) > sSc(1)
    fil_reg(end+1:sSc(2),1)     =   0;
end
Rc                              =   Vc*diag(fil_reg)*Vc';
dRc                             =   diag(Rc);
modv_Qc(modv_Qc(:,5)==1,10)     =   dRc; %diag of resolution matrix

% WE ALSO OUTPUT THE DIAGONAL OF THE RESOLUTION MATRIX
% using the filter functions for Q
sS                              =   size(S);
fil_reg                         =   fil_fac(diag(S),tik0_reg);
if sS(2) > sS(1)
    fil_reg(end+1:sS(2),1)      =   0;
end
R                               =   V*diag(fil_reg)*V';
dR                              =   diag(R);
modv_Q(modv_Q(:,5)==1,10)       =   dR; %diag of resolution matrix

%% Save in and outside Murat
Murat.data.modvPeakDelay        =   modv_pd;
Murat.data.modvQc               =   modv_Qc;
Murat.data.modvQ                =   modv_Q;

% save peak-delay
FName                           =   'peakdelay.txt';
save(fullfile(FPath, FLabel, FName), 'modv_pd','-ascii');

% save Qc
FName                           =   'Qc.txt';
save(fullfile(FPath, FLabel, FName), 'modv_Qc','-ascii');

% save Q
FName                           =   'Q.txt';
save(fullfile(FPath, FLabel, FName), 'modv_Q','-ascii');
