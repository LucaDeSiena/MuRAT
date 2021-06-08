% function [t0,A0,N,coda,t]=paasschens_function(r,v,B0,Le_1,dt,T)
% This fuction provides the Paasschens function for a fix r, with constants
% v,B0,Le_1, for points in the vector t until t_max given by T
% The Paasschens function is composed of a delta plus a coda
% The delta is in t0 = r/v
% The coda starts in t0 = r/v
% The amplitude of the delta is A0
% The coda is computed from t0 to T
% The coda is reported at samples starting at t0
% Only the N points for which coda is computed are provided
% Sample at t=t0 tends to infinity; since it is integrable, the
%   corresponding value is provided:  integral_0^T [E(r,t0+dt)dt] = E(r,t0+T) T 4/3; 
%   the size of the interval is half, so a value E(r,t0+dt)*2/3 is stored at t=t0 
%   This improves accuracy for not extremely small dt (by including sample
%   at t=t0).
% Input parameters:
%    r: distance to source
%    v: velocity
%    B0, Le_1: attenuation constants (scattering/intrincsic att. related constants) 
%    dt: time resolution
%    T: final time of interest
%
% Output parameters:
%    t0: time of the delta and beginning of coda: t0=r/v
%    A0: Amplitude of the delta
%    N: number of points from t0 to T (or from 0 to T-t0) for which coda is computed 
%    coda: the samples of the coda
%    t: time vector

function [t0,A0,N,coda,t]=paasschens_function_v1(r,v,B0,Le_1,dt,T)

t0=r/v;    % time of the delta in s
t=(t0:dt:(T+2*dt))';
N=length(t);

A0=exp(-Le_1*v*t0)./(4*pi*r.^2*v);

f1=(1-t0.*t0./(t.*t)).^(1/8);
f2=(3*B0*Le_1./(4*pi*v*t)).^(3/2);
f3=exp(-Le_1*v*t);
x=v*t*B0*Le_1.*(1-t0.*t0./(t.*t)).^(3/4);
f4=F_function(x);
coda=f1.*f2.*f3.*f4;
coda(1)=coda(2)*2/3;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F=F_function(x)
cond=x>1e-30;
F=zeros(size(x));
y=x(cond);
F(cond)=exp(y).*sqrt(1+2.026./(y+1e-30));
%F=exp(x).*sqrt(1+2.026./(x+1e-30));  % adding 1e-30 to avoid NaN
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%