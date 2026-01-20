function [d0, Q0, c]   =   Murat_lsqlinQmean(tCm,tWm,Q_k,cf_k,D,l,...
    time0_k,rapsp_k)
% [d0, Q0, c]   =   Murat_lsqlinQmean(tCm,tWm,Q_k,cf_k,D,l,...
%   time0_k,rapsp_k)
%
% INVERTS with minimum least squares to obtain average Q
%
% Input parameters:
%    tCm:           starting coda time
%    tWm:           coda window length
%    Q_k:           average coda attenuation
%    cf_k:          central frequency
%    D:             diffusion constant
%    l:             total ray length in km
%    time0_k:       travel time per frequency
%    rapsp_k:       energy ratio per frequency
%
% Output parameters:
%    d1:            data for the inversion - variations from average
%    const_Qc_k:    constant obtained using the average source-station Qc
%    constQmean_k:  contains geometrical spreading, Q, uncertainties
%    equationQ:     equation to be compared with data in test

%%
% Data creation for the true inversion, removing the parameters
% pre-calculated using the diffusion model

tLapse  =   tCm+tWm/2;
c       =   (1.5*log(4*pi*D*tLapse)...
            + l.^2/4/D./tLapse...
            + 2*pi*Q_k*cf_k.*tLapse)...
            /2/pi/cf_k;

data    =   log(rapsp_k.*l.^2)/2/pi/cf_k;

d0      =   data-c;

G       =   -time0_k;

% Storing inverted parameters
Q0      =   lsqlin(G,d0);
Q0(2)   =   (G'*G)^(-1)*G'*cov(d0)^-1*G*(G'*G)^(-1);

end