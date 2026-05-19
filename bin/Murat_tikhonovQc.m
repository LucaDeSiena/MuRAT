function [mtik0C,residualQc_k,LcQc]   =...
    Murat_tikhonovQc(outputLCurve,Gc,bQm,dampValue,x0,maxIt)
% function [mtik0C,residualQc_k,LcQc]   =...
%   Murat_tikhonovQc(outputLCurve,Gc,bQm,dampValue)  
%
% INVERTS with weighted tikhonov and creates L-curve and data for Qc.
%
% Input parameters:
%    outputLCurve:  flag to output the L curve
%    Gc:            weighted inversion matrix
%    bQm:           weighted inverse Qc values
%    dampValue:     damping parameter
%
% Output parameters:
%    mtik0C:        inversion parameter
%    residualQc_k:  residuals for Qc inversion
%    LcQc:          figure of L curve for the Qc method

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
% inputs you have: Gc (matrix), dampValue (scalar), bQm (vector)
n       =   size(Gc,2);
Afun    = @(x) Gc * x;
Atfun   = @(y) Gc' * y;
Lfun    = @(x) x;        % identity regularizer
Ltfun   = @(y) y;
lambda  = dampValue;
opts.maxit      = maxIt;
opts.tol        = 1e-7;
opts.verbose    = true;

x   =   Murat_tikhonov_nonneg_fista(Afun, Atfun, Lfun, Ltfun, bQm, n,...
    lambda, opts);

obj             =   norm((bQm-Gc*x).^2);
residualQc_k    =   obj/obj0;
mtik0C          =   x;

end