function cCoda      =   Murat_codaConstant(cf_k,l,Qc_k,tCm,te,D,v,tl)
% function cCoda    =   Murat_codaConstant(cf_k,l,Qc_k,tCm,te,D,v,tl)
%
% COMPUTES the coda constant per ray from the diffusion model
%   (Wu 1985, Wu and Aki, 1988; Sato et al., 2012).
%
% The data are d = log(E*r^2/Ec)/(2*pi*f), where Murat_body measures E as
%   the direct energy integrated over the body window tl and Ec as the
%   coda energy integrated over [tCm te] and rescaled to tl, i.e.
%   Ec = tl*mean(E_D). With the direct energy W*exp(-2*pi*f*t/Q)/(4*pi*v*r^2)
%   and the diffusive energy E_D = W*(4*pi*D*t)^-1.5*
%   exp(-r^2/(4*D*t) - 2*pi*f*Qc*t), d = -t/Q + cCoda with
%
%   cCoda = [1.5*log(4*pi*D) - log(mean(t^-1.5*exp(...))) ...
%           - log(4*pi*v*tl)]/(2*pi*f)
%
% Input parameters:
%    cf_k:          central frequency
%    l:             total ray length in km
%    Qc_k:          inverse Qc per ray
%    tCm:           start of the coda window
%    te:            end of the coda window
%    D:             diffusion constant (scalar, or column vector for a
%                   grid search with l, Qc_k, tCm and te as row vectors)
%    v:             average S-wave velocity
%    tl:            length of the body-wave window
%
% Output parameters:
%    cCoda:         coda constant, NaN where it is undefined

% Log of the mean of t^-1.5*exp(...) over the window (trapezoidal rule),
%   computed with log-sum-exp to avoid underflow at large r^2/(4*D*t)
nT      =   51;
w       =   ones(1,nT);
w([1 nT])   =   0.5;
logF    =   @(t) -1.5*log(t) - l.^2./(4*D.*t) - 2*pi*cf_k*Qc_k.*t;

Xmax    =   -Inf;
for j = 1:nT
    t   =   tCm + (te-tCm)*(j-1)/(nT-1);
    Xmax=   max(Xmax, logF(t));
end
S       =   0;
for j = 1:nT
    t   =   tCm + (te-tCm)*(j-1)/(nT-1);
    S   =   S + w(j)*exp(logF(t) - Xmax);
end
logA    =   Xmax + log(S) - log(nT-1);

cCoda   =   (1.5*log(4*pi*D) - logA - log(4*pi*v*tl))/2/pi/cf_k;
cCoda(~isfinite(cCoda)) =   NaN;
end
