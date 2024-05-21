function [d1, const_Qc_k, constQmean_k, equationQ]   =...
    Murat_lsqlinQmean(tCm,tWm,Q_k,cf_k,sped,luntot_k,time0_k,rapsp_k)
% function [d1, const_Qc_k, constQmean_k, equationQ]   =...
%     Murat_lsqlinQmean(tCm,tWm,Q_k,cf_k,sped,luntot_k,time0_k,rapsp_k)
%
% INVERTS with weighted tikhonov and creates L-curve and data for Q.
%
% Input parameters:
%    tCm:           lapse time
%    tWm:           coda window
%    Q_k:           inverse Qc
%    cf_k:          central frequency
%    sped:          spectral decay
%    luntot_k:      total ray length
%    time0_k:       total travel time
%    rapsp_k:       energy ratio
%    
% Output parameters:
%    d1:            data for average Q inversion
%    const_Qc_k:    constant for average Q inversion depending on Qc
%    constQmean_k:  average geometrical spreading and Q
%    equationQ:     sum of all three terms of the CN equation
%
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