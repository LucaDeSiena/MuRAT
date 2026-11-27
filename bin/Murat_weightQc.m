function w = Murat_weightQc(uncertaintyQc, QcMeasurement)
% function w = Murat_weightQc(uncertaintyQc, QcMeasurement)
%
% RETURNS the weight of each inverse Qc measurement, larger for better
% measurements. Murat_Qc stores a different quantity for each method:
%   'Linearized': the absolute correlation coefficient |R| of the fit,
%                 so w = |R|^2;
%   'NonLinear':  the uncertainty sigma of the grid-search fit,
%                 so w = 1/sigma^2.
% Measurements with zero, negative or NaN uncertainty get zero weight.
%
% Input parameters:
%    uncertaintyQc:  Murat.Qc.uncertaintyQc (any size)
%    QcMeasurement:  'Linearized' or 'NonLinear'
%
% Output parameters:
%    w:              weights, same size as uncertaintyQc
%
% Example:
%    w = Murat_weightQc([0.9; 0.5], 'Linearized');   % [0.81; 0.25]

u   =   uncertaintyQc;
ok  =   isfinite(u) & u > 0;
w   =   zeros(size(u));

switch QcMeasurement
    case 'Linearized'
        w(ok)   =   u(ok).^2;
    case 'NonLinear'
        w(ok)   =   1 ./ u(ok).^2;
    otherwise
        error('Murat_weightQc: unknown QcMeasurement ''%s''.', ...
            QcMeasurement);
end

end
