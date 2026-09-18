function  retain_Qm_i    =  Murat_retainQc(fT,Qm_i,RZZ_i,QcM)
% function  retain_Qm_i  =  Murat_retainQc(fT,Qm_i,RZZ_i,Ac_i,QcM)
%
% Creates all constraints for Qc inversion
%
% Input parameters:
%    fT:    treshold on uncertainty
%    Qm_i:  values of coda attenuation
%    RZZ_i: values of uncertainty
%    Ac_i:  coda attenuation matrix
%    QcM:   Linearized or Non Linear measurement
%
% Output parameters:
%    retain_Qm_i:       keeps tab on which waveforms are kept for imaging

if isequal(QcM,'Linearized')
    retainQmTemp    =   Qm_i>0 & RZZ_i<1/fT;
    retain_Qm_i     =   Qm_i>0 & RZZ_i<1/fT &...
        Qm_i < mean(Qm_i(retainQmTemp))+2*std(Qm_i(retainQmTemp));

elseif isequal(QcM,'NonLinear')
    retain_Qm_i     =   Qm_i>0 & RZZ_i<fT;
        
end

end