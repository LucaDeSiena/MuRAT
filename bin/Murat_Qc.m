function [inverseQc_i, uncertaintyQc_i]     =   Murat_Qc(cf,sped,sp_i,...
    cursorCodaStart_i,cursorCodaEnd_i,tCoda_i,srate_i,QcMeasurement)
% function [inverseQc_i, uncertaintyQc_i]	=   Murat_Qc(cf,sped,sp_i,...
%     cursorCodaStart_i,cursorCodaEnd_i,tCoda_i,srate_i,QcMeasurement)
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
%    QcMeasurement:     decides if Linearized or NonLinear solutions
% Output parameters:
%    inverseQc_i:       inverse coda attenuation factor
%    uncertaintyQc_i:   uncertainty on inverse coda attenuation factor

lcf                                         =   length(cf);
inverseQc_i                                 =   zeros(lcf,1);
uncertaintyQc_i                             =   zeros(lcf,1);

for i = 1:lcf
    
    envelopeC                               =...
        sp_i(cursorCodaStart_i:cursorCodaEnd_i,i);
    lEnvelopeC                              =   length(envelopeC);
    
    cf_i                                    =   cf(i);
    
    if isempty(envelopeC)
        
        inverseQc_i(i)                      =   0;
        uncertaintyQc_i(i)                  =   0;
        
    else
        
        if isequal(QcMeasurement,'Linearized')
            
            lapseT                          =...
                (tCoda_i+1/srate_i:1/srate_i:tCoda_i+lEnvelopeC/srate_i)';
            
            [linearFit, uncertaintyFit]     =...
                estimatesLinear(lapseT,envelopeC,lEnvelopeC,sped,cf_i);
            
            if linearFit(1)<0
                inverseQc_i(i)              =   -linearFit(1);
                uncertaintyQc_i(i)          =   abs(uncertaintyFit(1,2));
            else
                inverseQc_i(i)              =   0;
                uncertaintyQc_i(i)          =   0;
            end
            
        elseif isequal(QcMeasurement,'NonLinear')
            QValues                         =   0:10^-5:10^-1;
            
            lWindow                         =...
                round((cursorCodaEnd_i-cursorCodaStart_i)/srate_i);
            lapseT                          =...
                (tCoda_i+0.5:tCoda_i+lWindow-0.5)';
            
            [nonLinearFit, uncertaintyFit]  =   estimatesNonLinear(lapseT,...
                envelopeC,QValues,sped,lWindow,cf_i,srate_i);
            
            inverseQc_i(i)                  =   nonLinearFit;
            uncertaintyQc_i(i)              =   uncertaintyFit;
            
        else
            
            error('Unknown strategy to calculate Qc');
            
        end
        
    end
    
end

uncertaintyQc_i(inverseQc_i==0)             =   0;

end
%%
% Calculations in the linarized case.
function [linearFit, uncertaintyFit]        =...
    estimatesLinear(lapseT,envelopeC,lEnvelopeC,sped,cf_i)

%Only evaluate central time series
edgeC                                       =   floor(0.05*lEnvelopeC);
lapseTime                                   =   lapseT(edgeC:end-edgeC);
spcm1                                       =   envelopeC(edgeC:end-edgeC);

weigthEnergy                                =   spcm1.*lapseTime.^sped;
logWEnergy                                  =   log(weigthEnergy)/2/pi/cf_i;

linearFit                                   =...
    polyfit(lapseTime,logWEnergy,1);
uncertaintyFit                              =...
    corrcoef([lapseTime,logWEnergy]);

end

%%
%%
% Calculations in the non-linar (grid-search) case.

function [nonLinearFit, uncertaintyFit]     =...
    estimatesNonLinear(lapseT,envelopeC,QValues,sped,lWindow,cf_i,srate_i)

lQValues                                    =   length(QValues);

dObs                                        =   zeros(lWindow,1);

for k = 1:lWindow
    ntm                                     =   (k-1)*srate_i + 1:k*srate_i;
    dObs(k,1)                               =   mean(envelopeC(floor(ntm)));
end

dObserved                                   =   dObs(1:end-1)/dObs(end);

E                                           =   zeros(lQValues,1);
for n = 1:lQValues
    dPre                                    =...
        lapseT.^(-sped).*exp(-2*pi*cf_i.*lapseT*QValues(n));
    
    dPredicted                              =   dPre(1:end-1)/dPre(end);
    E(n,1)                                  =...
        sum(abs(dObserved-dPredicted));
end

[Emin, indexEmin]                           =   min(E);
nonLinearFit                                =   QValues(indexEmin);
uncertaintyFit                              =   1/Emin;

end

