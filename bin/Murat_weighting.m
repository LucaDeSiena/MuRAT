function Wc         =   Murat_weighting(RZZ_k,QcM)
%
% Input parameters:
%   RZZ:          uncertainty on Qc
%   QcM:          Linearized or Non Linear measurements.
% Output parameters:
%   Wc:           weighting matrix

if isequal(QcM,'Linearized')
    
    W1              =   RZZ_k<0.3;
    W2              =   RZZ_k<0.5;
    W3              =   RZZ_k<0.7;
    W4              =   W3+W2+W1;
    Wc              =   diag(1./(W4+1));
elseif isequal(QcM,'NonLinear')
    probDist        =   fitdist(RZZ_k,'LogNormal');
    confInterval    =   paramci(probDist);
    
    sortRZZ_k       =   sort(RZZ_k);
    pdfValue        =...
        pdf(probDist,sortRZZ_k)/max(pdf(probDist,sortRZZ_k));
    
    W1              =   pdfValue<0.3;
    W2              =   pdfValue<0.5;
    W3              =   pdfValue<0.7;
    W4              =   RZZ_k<confInterval(1) & RZZ_k>confInterval(2);
    W4              =   W4+W3+W2+W1;
    Wc              =   diag(1./(W4+1));

end

end