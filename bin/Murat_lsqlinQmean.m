function [d1, const_Qc_k, constQmean_k, equationQ]   =...
    Murat_lsqlinQmean(tCm,tWm,Q_k,cf_k,sped,luntot_k,time0_k,rapsp_k)
% function [d1, const_Qc_k, constQmean_k, equationQ]   =...
%     Murat_lsqlinQmean(tCm,tWm,Q_k,cf_k,sped,luntot_k,time0_k,rapsp_k)
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

%%
% Data creation for the true inversion, removing the pre-calculated
% parameters.
const_Qc_k                      =   (tCm+tWm/2).^-sped...
                                    .*exp(-2*pi*Q_k*cf_k.*(tCm+tWm/2));

d0                              =...
                    log(rapsp_k)/2/pi/cf_k + log(const_Qc_k)/2/pi/cf_k;

G                               =   -log(luntot_k)/2/pi/cf_k;
G(:,2)                          =   -time0_k;

cova                            =   (G'*G)^(-1)*G'*cov(d0)*G*(G'*G)^(-1);

constQmean_k(:,1)               =   lsqlin(G,d0(:,1));
constQmean_k(:,2)               =   sqrt(diag(cova));

d1                              =	d0  + ...
                            constQmean_k(1,1)*log(luntot_k)/2/pi/cf_k...
                            + time0_k*constQmean_k(2,1);

equationQ                       =   -log(const_Qc_k)...
                                    -constQmean_k(1,1)*log(luntot_k)...
                                    -2*pi*cf_k*time0_k*constQmean_k(2,1);

end