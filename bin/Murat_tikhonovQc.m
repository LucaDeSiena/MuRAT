function [mtik0C,residualQc_k,LcQc,tik0_regC]   =...
    Murat_tikhonovQc(outputLCurve,Gc,bQm,lCurveQc)
% function [mtik0C,residualQc_k,LcQc,tik0_regC]   =...
%     Murat_tikhonovQc(outputLCurve,Gc,bQm,lCurveQc)   
%
% INVERTS with weighted tikhonov and creates L-curve and data for Qc.
%
% Input parameters:
%    outputLCurve:  flag to output the L curve
%    Qm_k:          inverse Qc values
%    Wc:            weighting matrix
%    Gc:            inversion matrix
%    lCurveQc:      damping parameter input from start
%
% Output parameters:
%    mtik0C:        inversion parameter
%    residualQc_k:  residuals for Qc inversion
%    LcQc:          figure of L curve for the Qc method
%    tik0_regC:     regularization parameters

[Uc,Sc,Vc]                                      =   svd(Gc);

LcQc                                            =...
    figure('Name','L-curve Qc','NumberTitle','off');
[rho,eta,reg_param]                             =...
    l_curve_tikh_svd(Uc,diag(Sc),bQm,100);
plot_lc(rho,eta,'-',1,reg_param)

if outputLCurve == 1
    tik0_regC                                   =...
        input('Your damping parameter for coda: ');
else
    tik0_regC                                   =   lCurveQc;
end

mtik0C                                          =...
    tikhonov(Uc,diag(Sc),Vc,bQm,tik0_regC);
residualQc_k                                    =...
    sum(abs(bQm-Gc*mtik0C).^2);
end