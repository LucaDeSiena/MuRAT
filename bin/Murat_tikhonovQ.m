function [mtik0,residualQ_k,LcCN] =...
    Murat_tikhonovQ(outputLCurve,A,d1,dampValue,x0,maxIt)
% function [mtik0,residualQ_k,LcCN]=...
%     Murat_tikhonovQ(outputLCurve,A,d1,dampValue)
%
% INVERTS with weighted tikhonov and creates L-curve and data for Q.
%
% Input parameters:
%    outputLCurve:  flag to output the L curve
%    A:             CN inversion matrix
%    d1:            weighted data
%    dampValue:     damping parameter
%
% Output parameters:
%    mtik0:         inversion parameter
%    residualQ_k:   residuals for Q inversion
%    LcCN:          figure of L curve for the CN method

[U,S,~]         =   svd(A);

if outputLCurve == 1
    visibility  = 'on';
else
    visibility  = 'off';
end

LcCN = figure('Name','L-curve Q','NumberTitle','off','visible',visibility);

[rho,eta,reg_param]     =   l_curve_tikh_svd(U,diag(S),d1,100);
plot_lc(rho,eta,'-',1,reg_param)
xlabel('Regularization Parameter (rho)');
ylabel('Residual Norm (eta)');
grid on;

obj0    =   norm((d1-A*x0).^2);

n       = size(A,2);
Afun    = @(x) A * x;      % forward: y = A*x
Atfun   = @(y) A' * y;    % adjoint:  x = A'*y
Lfun    = @(x) x;        % identity regularizer
Ltfun   = @(y) y;
lambda  = dampValue;
opts.maxit      = maxIt;
opts.tol        = 1e-7;
opts.verbose    = true;

x   =   Murat_tikhonov_nonneg_fista(Afun, Atfun, Lfun, Ltfun, d1, n,...
    lambda, opts);

mtik0       =   x;

obj         =   norm((d1-A*mtik0).^2);
residualQ_k =   obj/obj0;

end