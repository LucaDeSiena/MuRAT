function sigma = Murat_uncertaintyNonLinear(QValues, Evec, nData)
% function sigma = Murat_uncertaintyNonLinear(QValues, Evec, nData)
%
% RETURNS the uncertainty of a 'NonLinear' (grid-search) inverse Qc
% measurement from the experimental probability density of the misfit.
%
% The grid search minimises the L1 misfit E(q) = sum(|dObs - dPre(q)|).
% An L1 misfit corresponds to a Laplace likelihood
%       p(q) ~ exp(-(E(q) - Emin)/b),
% whose scale b is estimated from the best fit as b = Emin/nData
% (maximum-likelihood estimate). sigma is the standard deviation of
% inverse Qc under p(q). It grows with the misfit of the fit and does not
% depend on the amplitude scale of the data. It is never smaller than the
% grid step, so a perfect fit still gets a finite weight.
%
% Input parameters:
%    QValues:   inverse Qc values tested by the grid search (uniform grid)
%    Evec:      L1 misfit for each value of QValues
%    nData:     number of data points entering the misfit
%
% Output parameters:
%    sigma:     uncertainty on inverse Qc

QValues     =   QValues(:);
Evec        =   Evec(:);
% Misfits that are NaN (predictions underflow to 0 for large trial values,
% late coda at high frequency) have zero probability
Evec(~isfinite(Evec))   =   Inf;
dq          =   min(diff(QValues));
Emin        =   min(Evec);

if ~isfinite(Emin) || nData < 1
    sigma   =   0;
    return
end

b           =   Emin/nData;
if b <= 0
    sigma   =   dq;
    return
end

p           =   exp(-(Evec - Emin)/b);
p           =   p/sum(p);
qMean       =   sum(p.*QValues);
sigma       =   sqrt(sum(p.*(QValues - qMean).^2));
sigma       =   max(sigma, dq);

end
