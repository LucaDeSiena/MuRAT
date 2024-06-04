function [d1, const_Qc_k, constQmean_k, equationQ]   =...
    Murat_lsqlinQmean(tCm,tWm,Q_k,cf_k,sped,luntot_k,time0_k,rapsp_k)
% function [d1, const_Qc_k, constQmean_k, equationQ]   =...
%     Murat_lsqlinQmean(tCm,tWm,Q_k,cf_k,sped,luntot_k,time0_k,rapsp_k)
%
% INVERTS with minimum least squares to obtain average Q
%
% Input parameters:
%    tCm:           starting coda time
%    tWm:           coda window length
%    Q_k:           average coda attenuation
%    cf_k:          central frequency
%    sped:          spectral decay
%    luntot_k:      total length per frequency
%    time0_k:       travel time per frequency
%    rapsp_k:       energy ratio per frequency
%
% Output parameters:
%    d1:            data for the inversion - variations from average
%    const_Qc_k:    constant obtained using the average source-station Qc
%    constQmean_k:  contains geometrical spreading, Q, uncertainties
%    equationQ:     equation to be compared with data in test

%%
% Data creation for the true inversion, removing the pre-calculated
% parameters.

% Include info on Qc
const_Qc_k                      =   (tCm+tWm/2).^-sped...
                                    .*exp(-2*pi*Q_k*cf_k.*(tCm+tWm/2));
% Data vector for inversion of average Q
d0                              =...
                    log(rapsp_k)/2/pi/cf_k + log(const_Qc_k)/2/pi/cf_k;

G                               =   -log(luntot_k)/2/pi/cf_k;
G(:,2)                          =   -time0_k;

cova                            =   (G'*G)^(-1)*G'*cov(d0)*G*(G'*G)^(-1);

% Storing inverted parameters
constQmean_k(:,1)               =   lsqlin(G,d0(:,1));
constQmean_k(:,2)               =   sqrt(diag(cova));

% Calculate data vector for average Q and geometrical spreading
d1                              =	d0  + ...
                            constQmean_k(1,1)*log(luntot_k)/2/pi/cf_k...
                            + time0_k*constQmean_k(2,1);

% Equation for average Q to fit data
equationQ                       =   -log(const_Qc_k)...
                                    -constQmean_k(1,1)*log(luntot_k)...
                                    -2*pi*cf_k*time0_k*constQmean_k(2,1);

end