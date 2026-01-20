function [mtik0C,residualQc_k,LcQc]   =...
    Murat_tikhonovQc(outputLCurve,Gc,bQm,dampValue,x0)
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

[Uc,Sc,~]          =   svd(Gc);

if outputLCurve == 1
    visibility      = 'on';
else
    visibility      = 'off';
end

LcQc = figure('Name','L-curve Qc','NumberTitle','off','visible',visibility);

[rho,eta,reg_param] =   l_curve_tikh_svd(Uc,diag(Sc),bQm,100);
plot_lc(rho,eta,'-',1,reg_param)
xlabel('Regularization Parameter (rho)');
ylabel('Residual Norm (eta)');
grid on;

obj0                =   norm((bQm-Gc*x0).^2);
Gaug                =   [Gc; dampValue*eye(numel(x0))];
baug                =   [bQm; dampValue*x0];
mtik0C              =   lsqnonneg(Gaug, baug);
obj                 =   norm((bQm-Gc*mtik0C).^2);

%mtik0C             =   tikhonov(Uc,diag(Sc),Vc,bQm,dampValue);
residualQc_k        =   obj/obj0;

end