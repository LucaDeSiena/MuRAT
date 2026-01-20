function [mtik0,residualQ_k,LcCN] =...
    Murat_tikhonovQ(outputLCurve,A,d1,dampValue,x0)
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

[rho,eta,reg_param] =   l_curve_tikh_svd(U,diag(S),d1,100);
plot_lc(rho,eta,'-',1,reg_param)
xlabel('Regularization Parameter (rho)');
ylabel('Residual Norm (eta)');
grid on;

obj0                =   norm((d1-A*x0).^2);
Gaug                =   [A; dampValue*eye(numel(x0))];
baug                =   [d1; dampValue*x0];
mtik0               =   lsqnonneg(Gaug, baug);
obj                 =   norm((d1-A*mtik0).^2);

%mtik0              =   tikhonov(U,diag(S),V,d1,dampValue);
residualQ_k         =   obj/obj0;

end