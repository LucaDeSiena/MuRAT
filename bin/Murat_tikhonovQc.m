function [mtik0C,residualQc_k,LcQc]   =...
    Murat_tikhonovQc(outputLCurve,Gc,bQm,dampValue,x0,maxIt,smoothValue,D)
% function [mtik0C,residualQc_k,LcQc]   =...
%   Murat_tikhonovQc(outputLCurve,Gc,bQm,dampValue,x0,maxIt,smoothValue,D)
%
% INVERTS with weighted tikhonov and creates L-curve and data for Qc.
%
% Minimises norm(Gc*m-bQm)^2 + dampValue^2*norm(m)^2 +
% smoothValue^2*norm(D*m)^2 with m >= 0. Without smoothValue and D (or
% with smoothValue = 0) only the damping term is used.
%
% Input parameters:
%    outputLCurve:  flag to output the L curve
%    Gc:            weighted inversion matrix
%    bQm:           weighted inverse Qc values
%    dampValue:     damping parameter
%    x0:            starting model (FISTA starts from it, negative
%                   values set to 0; also normalises the residual)
%    maxIt:         maximum FISTA iterations
%    smoothValue:   (optional) weight of the first-difference regulariser
%    D:             (optional) first-difference operator for the nodes in
%                   Gc (see Murat_gradientOperator)
%
% Output parameters:
%    mtik0C:        inversion parameter
%    residualQc_k:  residuals for Qc inversion
%    LcQc:          figure of L curve for the Qc method

if nargin < 7, smoothValue = []; end
if nargin < 8, D = []; end

[Uc,Sc,~]       =   svd(Gc);

if outputLCurve == 1
    visibility  = 'on';
else
    visibility  = 'off';
end

LcQc    =   figure('Name','L-curve Qc','NumberTitle','off','visible',visibility);

[rho,eta,reg_param]     =   l_curve_tikh_svd(Uc,diag(Sc),bQm,100);
plot_lc(rho,eta,'-',1,reg_param)
xlabel('Regularization Parameter (rho)');
ylabel('Residual Norm (eta)');
grid on;

obj0    =   norm((bQm-Gc*x0).^2);
n       =   size(Gc,2);
Afun    = @(x) Gc * x;
Atfun   = @(y) Gc' * y;
[Lfun, Ltfun, lambda] = Murat_regulariser(n, dampValue, smoothValue, D);
opts.maxit      = maxIt;
opts.tol        = 1e-7;
opts.verbose    = true;
opts.x0         = x0;

x   =   Murat_tikhonov_nonneg_fista(Afun, Atfun, Lfun, Ltfun, bQm, n,...
    lambda, opts);

obj             =   norm((bQm-Gc*x).^2);
residualQc_k    =   obj/obj0;
mtik0C          =   x;

end
