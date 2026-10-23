function [d0, Q0, mDest, cRay]   =   Murat_lsqlinQmean(cf_k,l,time0_k,...
    rapsp_k, Qc_k, tCm, te, v, tl, Dfix)
% [d0, Q0, mDest, cRay]   =   Murat_lsqlinQmean(cf_k,l,time0_k,rapsp_k,...
%    Qc_k, tCm, te, v, tl, Dfix)
%
% INVERTS with minimum least squares to obtain average Q
%
% The data d = log(E*r^2/Ec)/(2*pi*f) follow d = -t*Q^-1 + cCoda, where
%   the coda constant of each ray depends on the diffusion constant D
%   (Murat_codaConstant). For a given D, the average inverse Q is the
%   least-squares solution after removing cCoda.
%
% Input parameters:
%    cf_f:          central frequency
%    l:             total ray length in km
%    time0_k:       travel time per frequency
%    rapsp_k:       energy ratio per frequency
%    Qc_k:          inverse Qc per ray
%    tCm:           start of the coda window
%    te:            end of the coda window
%    v:             average S-wave velocity
%    tl:            length of the body-wave window
%    Dfix:          (optional) theoretical diffusion constant. If empty
%                   or absent, D is estimated as the value minimizing the
%                   misfit; otherwise D is fixed to Dfix and only the
%                   average inverse Q is inverted
%
% Output parameters:
%    d0:            data for the average inversion (coda constant removed)
%    Q0:            [mean coda constant, 0; average inverse Q, variance]
%    mDest:         estimated diffusion constant (Dfix if given)
%    cRay:          coda constant per ray

if nargin < 10
    Dfix    =   [];
end

data    =   log(rapsp_k.*l.^2)/2/pi/cf_k ;

% Misfit of the average inverse Q fit for a diffusion constant D
codaFit =   @(D) fitQ(data,time0_k,...
            Murat_codaConstant(cf_k,l,Qc_k,tCm,te,D,v,tl));

if ~isempty(Dfix)
    mDest   =   Dfix;
else
    % Grid search over log10(D), then local refinement
    logDs   =   linspace(-3, 3, 241);
    misfit  =   zeros(size(logDs));
    for i = 1:numel(logDs)
        misfit(i)   =   codaFit(10^logDs(i));
    end
    if all(isinf(misfit))
        error(['The coda constant is undefined for some rays for all '...
            'diffusion constants at %g Hz.'],cf_k)
    end
    [~, iMin]   =   min(misfit);
    if iMin == 1 || iMin == numel(logDs)
        warning(['The estimated diffusion constant at %g Hz is at the '...
            'edge of the search range [1e-3 1e3].'],cf_k)
    end
    iL      =   max(iMin-1, 1);
    iR      =   min(iMin+1, numel(logDs));
    logD    =   fminbnd(@(x) codaFit(10^x), logDs(iL), logDs(iR));
    mDest   =   10^logD;
end

cRay    =   Murat_codaConstant(cf_k,l,Qc_k,tCm,te,mDest,v,tl);
if any(isnan(cRay))
    error(['The coda constant is undefined for %d rays with the '...
        'diffusion constant %.4g at %g Hz.'],sum(isnan(cRay)),mDest,cf_k)
end

[~, QiEst, varQi]   =   fitQ(data,time0_k,cRay);

Q0      =   [mean(cRay) 0; QiEst varQi];
d0      =   data - cRay;
end

function [misfit, QiEst, varQi]  =   fitQ(data,time0_k,cRay)
% Least-squares average inverse Q for d - cRay = -t*Q^-1
if any(isnan(cRay))
    misfit  =   Inf;
    QiEst   =   NaN;
    varQi   =   NaN;
    return
end
dataQ   =   data - cRay;
G       =   -time0_k;
QiEst   =   (G'*dataQ)/(G'*G);
res     =   G*QiEst - dataQ;
misfit  =   res'*res;
dof     =   max(numel(G) - 1, 1);
varQi   =   misfit/dof/(G'*G);
end
