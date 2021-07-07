function [sp_i,SAChdr_i,srate_i]...
                        =   Murat_envelope(Murat,listSac_i)
% ENVELOPE calculation, the code also outputs the picking time 

cf                      =   Murat.input.centralFrequency;
lcf                     =   length(cf);
% Imports SAC files
[~,sisma,SAChdr]        =   fget_sac(listSac_i); 
% Sampling rate
srate_i                 =   1/SAChdr.times.delta;
% Remove mean and trends
sisma                   =   detrend(sisma,1);
% Seismogram length
lsis                    =   length(sisma);
% Tapering
tu                      =   tukeywin(lsis,0.05);
tsisma                  =   tu.*sisma;
sp_i                =   zeros(lsis,lcf);

%% Calculations
for i = 1:lcf
    % Filter creation - in loop for different frequencies
    Wn              =...
        ([cf(i)-cf(i)/3 cf(i)+cf(i)/3]/srate_i*2);
    [z,p,k]         =   butter(4,Wn,'bandpass');
    [sos,g]         =   zp2sos(z,p,k);
    
    fsisma          =   filtfilt(sos,g,tsisma);
    
    SAChdr_i        =   SAChdr;
    [sp,~]          =   envelope(fsisma,round(srate_i),'rms');
    sp_i(:,i)       =   sp.^2;
end
end
