function [mtik0,residualQ_k,LcCN,tik0_reg]...
                                =...
    Murat_tikhonovQ(outputLCurve,A,d1,lCurveQ_k,flagShow)
% function [mtik0,residualQ_k,LcCN,tik0_reg,d1k,constQmean_k]...
%                                 =...
%     Murat_tikhonovQ(cfk,rtQ,outputLCurve,rapsp_k,const_Qc_k,...
%     luntot,time0,A,lCurveQ,flagShow)
%
% INVERTS with weighted tikhonov and creates L-curve and data for Q.
%
% Input parameters:
%    cfk:           central frequency
%    rtQ:           selected Q waveforms
%    outputLCurve:  flag to output the L curve
%    rapsp_k:       energy ratios
%    const_Qc_k:    constant from the average Qc
%    luntot:        total length
%    time0:         travel time
%    A:             CN inversion matrix
%    lCurveQ_k:     damping parameter input from start
%    flagShow:      flag to decide if show L-curve
%
% Output parameters:
%    mtik0:         inversion parameter
%    residualQ_k:   residuals for Q inversion
%    LcCN:          figure of L curve for the CN method
%    tik0_reg:      regularization parameters
%    d1:            data of the inversion
%    constQmean_k:  constant for the method

%%
% Creates L-curve figure as before and sets damping
[U,S,V]                         =   svd(A);
s                               =   diag(S);

if flagShow == 1
    LcCN                            =...
        figure('Name','L-curve Q','NumberTitle','off');
    [rho,eta,reg_param]             =...
        l_curve_tikh_svd(U,s,d1,100);
    plot_lc(rho,eta,'-',1,reg_param)
    
    if outputLCurve == 1
        tik0_reg                =...
            input('Your damping parameter for Q: ');
    else
        tik0_reg                =   lCurveQ_k;
    end
else
    tik0_reg                    =   lCurveQ_k;
    LcCN                        =   0;
end
mtik0                           =   tikhonov(U,s,V,d1,tik0_reg);
residualQ_k                     =   sum(abs(d1-A*mtik0).^2);

end