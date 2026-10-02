function [tempis,fsisma_i,SAChdr_i,srate_i]   =   Murat_sisma(cf,listSac_i)
% function [tempis,fsisma_i,SAChdr_i,srate_i]   =   Murat_sisma(cf,listSac_i)
%
% SEISMOGRAM calculation at all frequencies, outputs properties of seismograms 
%
% Input parameters:
%    cf:            central frequency
%    listSac_i:     list of SAC files
%    
% Output parameters:
%    tempis:        time vector of seismogram
%    fsisma_i:      seismogram filtered at different frequencies
%    SAChdr_i:      SAC header
%    srate_i:       sampling rate
 
lcf                     =   length(cf);
[tempis,sisma,SAChdr_i] =   fget_sac(listSac_i); 
srate_i                 =   1/SAChdr_i.times.delta;

sisma                   =   detrend(sisma,1);
lsis                    =   length(sisma);
tu                      =   tukeywin(lsis,0.05);
tsisma                  =   tu.*sisma;

fsisma_i                =   zeros(lsis,lcf);

if isequal(srate_i,-12345)
    error(['Waveform ' listaSac_i 'has no sampling rate!'])
end

for i = 1:lcf
    % Filter creation - in loop for different frequencies
    Wn                  =   ([cf(i)-cf(i)/3 cf(i)+cf(i)/3]/srate_i*2);
    [z,p,k]             =   butter(4,Wn,'bandpass');
    [sos,g]             =   zp2sos(z,p,k);
    
    fsisma_i(:,i)       =   filtfilt(sos,g,tsisma);
    
end
end
