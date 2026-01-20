function [inverseQc_i, uncertaintyQc_i] =   Murat_Qc(cf,...
  sp_i,cursorCodaStart_i,cursorCodaEnd_i,tCoda_i,srate_i,QcMeasurement)
% function [inverseQc_i, uncertaintyQc_i]   =   Murat_Qc(cf,...
% sp_i,cursorCodaStart_i,cursorCodaEnd_i,tCoda_i,srate_i,QcMeasurement)
%
% MEASURES Qc and its uncertainty.
%
% Input parameters:
%    cf:                central frequency
%    sp_i:              filtered envelope
%    cursorCodaStart_i: start of the coda on trace
%    cursorEndStart_i:  end of the coda on trace
%    tCoda_i:           start of the coda in seconds
%    srate_i:           sampling rate
%    QcMeasurement:     decides if Linearized or NonLinear solutions
%
% Output parameters:
%    inverseQc_i:       inverse coda attenuation factor
%    uncertaintyQc_i:   uncertainty on inverse coda attenuation factor

lcf                     =   numel(cf);
inverseQc_i             =   zeros(lcf,1);
uncertaintyQc_i         =   zeros(lcf,1);
srate_i                 =   round(srate_i);
si                      =   1/srate_i;

for i = 1:lcf
    envelopeC           =   sp_i(cursorCodaStart_i:cursorCodaEnd_i,i);
    lEnvelopeC          =   numel(envelopeC);
    cf_i                =   cf(i);

    switch QcMeasurement

        case 'Linearized'
            
            lapseT          =  (tCoda_i+si:si:tCoda_i+lEnvelopeC*si)';

            edgeC           =   floor(0.05*lEnvelopeC);
            lapseTime       =   lapseT(edgeC:end-edgeC);
            spcm1           =   envelopeC(edgeC:end-edgeC);

            weigthEnergy    =   spcm1.*lapseTime.^1.5;
            logWEnergy      =   log(weigthEnergy)/2/pi/cf_i;

            p               =   polyfit(lapseTime,logWEnergy,1);
            R               =   corrcoef([lapseTime,logWEnergy]);

            if p(1)<0
                inverseQc_i(i)      =   -p(1);
                RZZ_k       =   abs(R(1,2));
                w           =   1 ./ RZZ_k;
                w           = (w - min(w(:))) ./ (max(w(:)) - min(w(:)));
                uncertaintyQc_i(i)  =   w;
                
            else
                inverseQc_i(i)      =   0;
                uncertaintyQc_i(i)  =   0;
            end

        case 'NonLinear'
            QValues     =   10^-5:10^-5:10^-1;
            lWindow     =   round((cursorCodaEnd_i-cursorCodaStart_i)*si);

            % returns srate_i-by-nCols matrix
            src = envelopeC(:);
            nCols = min(lWindow, floor(numel(src) / srate_i));% n full col.
            if nCols == 0
                blk = [];   % no full blocks available
            else
                blk = reshape(src(1 : nCols * srate_i), srate_i, nCols);
            end
            
            dObs        =   mean(blk)';

            % normalized observed (exclude last sample, normalize by last)
            dObserved   =     dObs(1:end-1) / dObs(end);

            % lapse times for blocks
            lapseT_blk  =   (tCoda_i + 0.5 : tCoda_i + size(blk,2) - 0.5)';

            % prepare predicted model matrix for all QValues at once
            L           =   lapseT_blk(:);           
            const       =   2 * pi * cf_i; 

            % compute matrix of exponents: (lapseT * QValues') is m x nQ
            expMat      =   exp(- const * (L * QValues(:)'));  
            dPreMat     =   (L .^ -1.5) .* expMat;

            % normalize each column by its last element, then drop last row
            den         =   dPreMat(end, :);
            dPredicted  =   dPreMat(1:end-1, :) ./ den;

            % compute L1 misfit for each QValue
            Evec        =   sum(abs(dObserved - dPredicted), 1);

            % pick best fit and compute uncertainty
            [Emin, idx]     =   min(Evec);
            nonLinearFit    =   QValues(idx);
            uncertaintyFit  =   -log(Emin / mean(lapseT_blk)^1.5) /...
                (2 * pi * cf_i * mean(lapseT_blk));

            inverseQc_i(i)              =   nonLinearFit;
            uncertaintyQc_i(i)          =   uncertaintyFit;

        otherwise

            error('Unknown strategy to calculate Qc');
    end

end

% force zero uncertainty where inverseQc is zero
uncertaintyQc_i(inverseQc_i <= 0)       =     0;
inverseQc_i(uncertaintyQc_i<0)          =    0;

end

