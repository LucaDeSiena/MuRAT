function [mtik0,residualQ_k,LcCN,tik0_reg,d1k,constQmean_k]...
                                =...
    Murat_tikhonovQ(cfk,rtQ,outputLCurve,rapsp_k,const_Qc_k,...
    luntot,time0,A,lCurveQ,flagShow)
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
%    lCurveQ:       damping parameter input from start
%    flagShow:      flag to decide if show L-curve
%
% Output parameters:
%    mtik0:         inversion parameter
%    residualQ_k:   residuals for Q inversion
%    LcCN:          figure of L curve for the CN method
%    tik0_reg:      regularization parameters
%    d1k:           data of the inversion
%    constQmean_k:  constant for the method

d0                              =...
    log(rapsp_k)/2/pi/cfk + log(const_Qc_k)/2/pi/cfk;
G                               =   -log(luntot(rtQ))/2/pi/cfk;
G(:,2)                          =   -time0(rtQ);

%%
% The three parameters are the constant, the geometrical spreading, and
% the average Q and they are contained in constQmean. Computed with
% a least square inversion. Using the covariance matrix for
% uncertainties.
constQmean_k(:,1)               =   lsqlin(G,d0(:,1));
cova                            =...
    (G'*G)^(-1)*G'*cov(d0)*G*(G'*G)^(-1);
constQmean_k(:,2)               =   sqrt(diag(cova));


%%
% Data creation for the true inversion, removing the pre-calculated
% parameters.
d1k                             =...
    d0  + constQmean_k(1,1)*log(luntot(rtQ))/2/pi/cfk...
    + time0(rtQ)*constQmean_k(2,1);
%%
% Same procedure for inversion for Q.
[U,S,V]                         =   svd(A);

%%
% Creates L-curve figure as before and sets damping
if flagShow == 1
    LcCN                            =...
        figure('Name','L-curve Q','NumberTitle','off');
    [rho,eta,reg_param]             =...
        l_curve_tikh_svd(U,diag(S),d1k,100);
    plot_lc(rho,eta,'-',1,reg_param)
    
    if outputLCurve == 1
        tik0_reg                =...
            input('Your damping parameter for Q: ');
    else
        tik0_reg                =   lCurveQ;
    end
else
   tik0_reg                     =   lCurveQ;
   LcCN                         =   0;
end
mtik0                           =   tikhonov(U,diag(S),V,d1k,tik0_reg);
residualQ_k                     =   sum(abs(d1k-A*mtik0).^2);

end