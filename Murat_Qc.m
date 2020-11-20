function [inverseQc_i,uncertaintyQc_i]    =...
    Murat_Qc(Murat,sp_i,cursorCodaStart_i,cursorCodaEnd_i,tCoda_i,srate_i)

sped                    =   Murat.input.spectralDecay;
cf                      =   Murat.input.centralFrequency;
tWm                     =   Murat.input.codaWindow;

%Use shorter window if waveform is cut
if cursorCodaEnd_i > length(sp_i)
    spcm                =   sp_i(cursorCodaStart_i:length(sp_i));
else
    spcm                =   sp_i(cursorCodaStart_i:cursorCodaEnd_i);
end
%Compute time vector
tm                      =...
    (tCoda_i+1/srate_i:1/srate_i:tCoda_i+length(spcm)/srate_i)';

%Only evaluate central time series
edgeno                  =   floor(0.05*length(spcm));
tm1                     =   tm(edgeno:end-edgeno);
spcm1                   =   spcm(edgeno:end-edgeno);

%Using the linear approximation
if Murat.input.nonLinear == 0
    
    %THIS IS THE DATA VECTOR WITH LINEARISED THEORY
    
    EWz                 =   spcm1.*tm1.^sped; %Calculating energy
    lspmz               =   log(EWz)/2/pi/cf; % source-station data
    Rz                  =   corrcoef([tm1,lspmz]); %sets uncertainty
    polyz               =   polyfit(tm1,lspmz,1); %Qc^-1
    
    if polyz(1)<0 %only where you have a positive Qc
        inverseQc_i     =   -polyz(1);
        uncertaintyQc_i =   abs(Rz(1,2));
    else
        inverseQc_i     =   0;
        uncertaintyQc_i =   0;
    end
    
elseif Murat.input.nonLinear == 1
    
    %Solving with the grid-search
    nW                  =   Murat.input.fitLengthWindows;
    ntW                 =   Murat.input.fitNumberWindows;
    L1                  =   Murat.input.fitNumber;
    m1a                 =   Murat.input.fitTrialQc;

    %THIS IS THE SYSTEM OF EQUATIONS FOR THE NON-LINEAR SOLUTION
    %USING THE LINEAR INVERSE QC AS STARTING MODEL
    
    tlapse              =   (tCoda_i+nW/2:nW:tCoda_i+tWm-nW/2)';
    d_obs               =   zeros(ntW,1);%set up the data vector
    
    %For all the smaller windows in play...
    for k = 1:ntW
        %take the samples in the window...
        ntm             =   (k-1)*nW*srate_i + 1:k*nW*srate_i;
        %and compute the average energy.
        d_obs(k,1)      =   mean(spcm(floor(ntm)));
    end
    
    %normalize by the last window
    d_obs1              =   d_obs(1:end-1)/d_obs(end);
    
    %PREDICTED BY THEORY
    E                   =   zeros(L1,1);
    for n = 1:L1
        %data predicted
        d_pre           =   tlapse.^(-sped).*exp(-2*pi*cf.*tlapse*m1a(n));
        d_pre1          =   d_pre(1:end-1)/d_pre(end);%normalized
        E(n,1)          =   sum(abs(d_obs1-d_pre1));%L1 residual
    end
    [Emin, indexE]      =   min(E);%find the minimum to set the right Qc
    inverseQc_i         =   m1a(indexE);%Qc set
    
    % Uncertainty from the minimum residual, the bigger it is the more
    % uncertain the model becomes
    uncertaintyQc_i     =   1/Emin;
end
