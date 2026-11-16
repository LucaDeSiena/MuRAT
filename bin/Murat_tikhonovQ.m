function [mtik0,residualQ_k,LcCN] =...
    Murat_tikhonovQ(outputLCurve,A,d1,dampValue,x0,maxIt,smoothValue,D)
% function [mtik0,residualQ_k,LcCN]=...
%     Murat_tikhonovQ(outputLCurve,A,d1,dampValue,x0,maxIt,smoothValue,D)
%
% INVERTS with weighted tikhonov and creates L-curve and data for Q.
%
% Minimises norm(A*m-d1)^2 + dampValue^2*norm(m)^2 +
% smoothValue^2*norm(D*m)^2 with m >= 0. Without smoothValue and D (or
% with smoothValue = 0) only the damping term is used.
%
% Input parameters:
%    outputLCurve:  flag to output the L curve
%    A:             CN inversion matrix
%    d1:            weighted data
%    dampValue:     damping parameter
%    x0:            starting model (FISTA starts from it, negative
%                   values set to 0; also normalises the residual)
%    maxIt:         maximum FISTA iterations
%    smoothValue:   (optional) weight of the first-difference regulariser
%    D:             (optional) first-difference operator for the nodes in
%                   A (see Murat_gradientOperator)
%
% Output parameters:
%    mtik0:         inversion parameter
%    residualQ_k:   residuals for Q inversion
%    LcCN:          figure of L curve for the CN method

if nargin < 7, smoothValue = []; end
if nargin < 8, D = []; end

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
[Lfun, Ltfun, lambda] = Murat_regulariser(n, dampValue, smoothValue, D);
opts.maxit      = maxIt;
opts.tol        = 1e-7;
opts.verbose    = true;
opts.x0         = x0;

x   =   Murat_tikhonov_nonneg_fista(Afun, Atfun, Lfun, Ltfun, d1, n,...
    lambda, opts);

mtik0       =   x;

obj         =   norm((d1-A*mtik0).^2);
residualQ_k =   obj/obj0;

end
