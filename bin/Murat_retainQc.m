function  [retain_Qm_i,ray_crosses_Qc_i]    =...
    Murat_retainQc(fT,Qm_i,RZZ_i,Ac_i,QcM)
% function  [retain_Qm_i,ray_crosses_Qc_i,QcM]    =...
%     Murat_retainQc(fT,Qm_i,RZZ_i,Ac_i)
%
% Creates all constraints for Qc inversion
%
% Input parameters:
%    fT:                treshold on uncertainty
%    Qm_i:              values of coda attenuation
%    RZZ_i:             values of uncertainty
%    Apd_i:             coda attenuation matrix
%    QcM:               Linearized or Non Linear measurement
%
% Output parameters:
%    retain_Qm_i:       keeps tab on which waveforms are kept for imaging
%    ray_crosses_Qc_i:  keeps tab on which rays are kept for imaging

if isequal(QcM,'Linearized')
    retainQmTemp                            =   Qm_i>0 & RZZ_i>fT;
    retain_Qm_i                             =   Qm_i>0 & RZZ_i>fT &...
        Qm_i < mean(Qm_i(retainQmTemp))+2*std(Qm_i(retainQmTemp));
    
    Ac_retain_Qc_i                          =   Ac_i(retain_Qm_i,:);
    Ac_retain_Qc_i(Ac_retain_Qc_i<10^(-4))  =   0;
    s_Qc                                    =   sum(Ac_retain_Qc_i);
    ray_crosses_Qc_i                        =   s_Qc > 0.1;

elseif isequal(QcM,'NonLinear')
    retain_Qm_i                             =   Qm_i>0;
    Ac_retain_Qc_i                          =   Ac_i(retain_Qm_i,:);
    Ac_retain_Qc_i(Ac_retain_Qc_i<10^(-4))  =   0;
    s_Qc                                    =   sum(Ac_retain_Qc_i);
    ray_crosses_Qc_i                        =   s_Qc > 0.1;
    
    
end

end