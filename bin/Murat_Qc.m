function [inverseQc_i, uncertaintyQc_i]     =   Murat_Qc(cf,sped,sp_i,...
    cursorCodaStart_i,cursorCodaEnd_i,tCoda_i,srate_i)
% function [inverseQc_i, uncertaintyQc_i]     =   Murat_Qc(cf,sped,sp_i,...
%     cursorCodaStart_i,cursorCodaEnd_i,tCoda_i,srate_i)
%
% MEASURES Qc and its uncertainty.
%
% Input parameters:
%    cf:                central frequency
%    sped:              spectral decay
%    sp_i:              envelopes
%    cursorCodaStart_i: start of the coda on trace
%    cursorEndStart_i:  end of the coda on trace
%    tCoda_i:           start of the coda in seconds
%    srate_i:           sampling rate
%
% Output parameters:
%    inverseQc_i:       inverse coda attenuation factor
%    uncertaintyQc_i:   uncertainty on inverse coda attenuation factor

lcf                                         =   length(cf);
inverseQc_i                                 =   zeros(lcf,1);
uncertaintyQc_i                             =   zeros(lcf,1);

for i = 1:lcf
    
    spcm                                    =...
        sp_i(cursorCodaStart_i:cursorCodaEnd_i,i);
    
    tm                                      =...
        (tCoda_i+1/srate_i:1/srate_i:tCoda_i+length(spcm)/srate_i)';
    
    %Only evaluate central time series
    edgeno                                  =   floor(0.05*length(spcm));
    tm1                                     =   tm(edgeno:end-edgeno);
    spcm1                                   =   spcm(edgeno:end-edgeno);
    
    EWz                                     =   spcm1.*tm1.^sped;
    lspmz                                   =   log(EWz)/2/pi/cf(i);
    Rz                                      =   corrcoef([tm1,lspmz]);
    polyz                                   =   polyfit(tm1,lspmz,1);
    
    if polyz(1)<0
        inverseQc_i(i)                      =   -polyz(1);
        uncertaintyQc_i(i)                  =   abs(Rz(1,2));
    else
        inverseQc_i(i)                      =   0;
        uncertaintyQc_i(i)                  =   0;
    end        
    
end