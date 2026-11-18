function Murat_printRegularisation(name,normA,dampValue,smoothTikh,smoothValue)
% function Murat_printRegularisation(name,normA,dampValue,smoothTikh,smoothValue)
%
% PRINTS the norm of the inversion matrix and the absolute weights that
%   the relative damping and Tikhonov smoothing values correspond to.
%   Values from input files written before damping became relative must
%   be divided by norm(A) to give the same inversion.
%
% Input parameters:
%    name:          'Q' or 'Qc'
%    normA:         largest singular value of the (weighted) matrix
%    dampValue:     relative damping
%    smoothTikh:    true if first-difference Tikhonov smoothing is used
%    smoothValue:   relative smoothing (used only if smoothTikh)

msg = sprintf(['%s inversion: norm(A) = %.4g, damping %.4g ' ...
    '(absolute %.4g)'], name, normA, dampValue, dampValue*normA);
if smoothTikh
    msg = [msg sprintf(', smoothing %.4g (absolute %.4g)', ...
        smoothValue, smoothValue*normA)];
end
fprintf('%s\n', msg);
end
