function [t0,A0,N,coda,t]   =   Murat_paasschensFunction(r,v,B0,Le_1,dt,T)
% function [t0,A0,N,coda,t] =   Murat_paasschensFunction(r,v,B0,Le_1,dt,T)
%
% CREATES the Paasschens function for a fixed r, with constants v,B0,Le_1,
%	for points in the vector t until t_max given by T.
%
% Structure:
% The Paasschens function is composed of a delta plus a coda:
% The delta is in t0 = r/v;
% The coda starts in t0 = r/v;
% The amplitude of the delta is A0;
% The coda is computed from t0 to T;
% The coda is reported at samples starting at t0;
% Only the N points for which coda is computed are provided.
%
% Sample at t=t0 tends to infinity. As it is integrable, the
%   corresponding value is provided:
%   
%   integral_0^T [E(r,t0+dt)dt] = E(r,t0+T) T 4/3; 
%   
%   the size of the interval is half, so E(r,t0+dt)*2/3 is stored at t=t0. 
%   This improves accuracy for small dt that are not too small by including
%   the sample at t=t0.
%
% Input parameters:
%    r:         distance to source
%    v:         velocity
%    B0:        albedo
%    Le_1:      inverse extinction length 
%    dt:        time resolution
%    T:         final time of interest
%
% Output parameters:
%    t0:        time of the delta and beginning of coda: t0=r/v
%    A0:        Amplitude of the delta
%    N:         number of points from t0 to T (or from 0 to T-t0) of coda
%    coda:      the samples of the coda
%    t:         time vector
% 
% Authors: De La Torre & Del Pezzo, used first in Del Pezzo et al. 2018,
% Geosciences.

t0          =   r/v;
t           =   (t0:dt:(T+2*dt))';
N           =   numel(t);

A0          =   exp(-Le_1*v*t0)./(4*pi*r.^2*v);

% Precompute recurring terms
t_inv       =   1 ./ t;
t2          =   t .* t;
t0_over_t   =   t0 .* t0 ./ t2;
one_minus   =   1 - t0_over_t;
one_minus(one_minus < 0) = 0;

f1          =   one_minus .^ (1/8);
const_f2    =   (3 * B0 .* Le_1) / (4 * pi * v);

f2          =   (const_f2 .* t_inv) .^ (3/2);

f3          =   exp(-Le_1*v.*t);

f1_pow      =   one_minus .^ (3/4);
x           =   v*t.*B0.*Le_1.*f1_pow;

% Inline F_function: F = exp(x) .* sqrt(1 + 2.026 ./ x) for x > tiny
F           =   zeros(size(x));
cond        =   x > 1e-30;
y           =   x(cond);
F(cond)     =   exp(y) .* sqrt(1 + 2.026 ./ y);

coda        =   f1.*f2.*f3.*F;
coda(1)     =   coda(2)*2/3;

end