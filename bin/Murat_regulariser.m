function [Lfun, Ltfun, lambda] = Murat_regulariser(n, dampValue,...
    smoothValue, D)
% function [Lfun, Ltfun, lambda] = Murat_regulariser(n, dampValue,...
%    smoothValue, D)
%
% BUILDS the regularisation operator used by Murat_tikhonov_nonneg_fista,
% which adds lambda^2*norm(L*m)^2 to the misfit.
%
% Without smoothing (smoothValue empty or 0, or D empty) L is the identity
% and lambda = dampValue, i.e. zero-order Tikhonov as in MuRAT 4.0.0.
% With smoothing L = [dampValue*I; smoothValue*D] and lambda = 1, so the
% penalty is dampValue^2*norm(m)^2 + smoothValue^2*norm(D*m)^2.
%
% Input parameters:
%    n:             number of model parameters
%    dampValue:     damping weight
%    smoothValue:   weight of the first-difference term
%    D:             first-difference operator (rows-by-n)
%
% Output parameters:
%    Lfun, Ltfun:   handles for L and L'
%    lambda:        scalar weight passed to the solver
%
% Example:
%    [Lfun, Ltfun, lambda] = Murat_regulariser(size(A,2), 0, 0.1, D);

if isempty(smoothValue) || smoothValue == 0 || isempty(D)
    Lfun    =   @(x) x;
    Ltfun   =   @(y) y;
    lambda  =   dampValue;
    return
end

if size(D,2) ~= n
    error('Murat_regulariser: D has %d columns, expected %d.', size(D,2), n);
end

nD          =   size(D,1);
dv          =   abs(dampValue);
sv          =   abs(smoothValue);
Lfun        =   @(x) [dv * x; sv * (D * x)];
Ltfun       =   @(y) dv * y(1:n) + sv * (D' * y(n+1:n+nD));
lambda      =   1;

end
