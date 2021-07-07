function  [retain_Qm_i,ray_crosses_Qc_i]    =...
    Murat_retainQc(nonlinear,fT,Qm_i,RZZ_i,Ac_i)
% Creates all constraints for Qc inversion

switch nonlinear
    
    case 0
        retainQmTemp                        =   Qm_i>0 & RZZ_i>fT;
        retain_Qm_i                         =...
            Qm_i>0 & RZZ_i>fT &...
            Qm_i < mean(Qm_i(retainQmTemp))+2*std(Qm_i(retainQmTemp));
        
    case 1
        retainQmTemp                        =   Qm_i>0 & RZZ_i<fT;
        retain_Qm_i                         =...
            Qm_i>0 & RZZ_i<fT &...
            Qm_i < mean(Qm_i(retainQmTemp))+2*std(Qm_i(retainQmTemp));
        
end

Ac_retain_Qc_i                              =   Ac_i(retain_Qm_i,:);
Ac_retain_Qc_i(Ac_retain_Qc_i<10^(-3))      =   0;
s_Qc                                        =   sum(Ac_retain_Qc_i);
ray_crosses_Qc_i                            =   s_Qc > 0.1;

end