function  retain_Qm_i    =  Murat_retainQc(fT,Qm_i,RZZ_i,QcM)
% function  retain_Qm_i  =  Murat_retainQc(fT,Qm_i,RZZ_i,QcM)
%
% Creates all constraints for Qc inversion
%
% Input parameters:
%    fT:    treshold on uncertainty
%    Qm_i:  values of coda attenuation
%    RZZ_i: |R| of the fit ('Linearized') or uncertainty ('NonLinear')
%    QcM:   Linearized or Non Linear measurement
%
% Output parameters:
%    retain_Qm_i:       keeps tab on which waveforms are kept for imaging

if isequal(QcM,'Linearized')
    % Murat_Qc stores |R|: keep fits with |R| above the threshold
    retainQmTemp    =   Qm_i>0 & RZZ_i>fT;
    retain_Qm_i     =   Qm_i>0 & RZZ_i>fT &...
        Qm_i < mean(Qm_i(retainQmTemp))+2*std(Qm_i(retainQmTemp));

elseif isequal(QcM,'NonLinear')
    switch fT
        case 0
            retain_Qm_i = Qm_i > 0 & RZZ_i < max(RZZ_i);
        otherwise
            retain_Qm_i = Qm_i > 0 & RZZ_i < fT;
    end
        
end

end