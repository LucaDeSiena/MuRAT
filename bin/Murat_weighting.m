function Wc             =   Murat_weighting(RZZ_k)
% CREATES weighting for the coda method
%
% Input parameters:
%   RZZ:          uncertainty on Qc
% Output parameters:
%   Wc:           weighting matrix

W1              =   RZZ_k<0.3;
W2              =   RZZ_k<0.5;
W3              =   RZZ_k<0.7;
W4              =   W3+W2+W1;
Wc              =   diag(1./(W4+1));

end