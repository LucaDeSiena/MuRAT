function [image, minimizeVector,infoVector,minimizeValue]   =...
    Murat_minimiseHybrid(outputLCurve,G,bQ,lCurveQc_k,name)
% function [image, minimizeVector,infoVector,minimizeValue]   =...
%     Murat_minimiseIR(G,bQ,name)
% MINIMIZES weighted least square solution with an hybrid iterative method
% based on generalised cross validation
% CALCULATES the damping value and PLOTS the costs.
%
% Input parameters:
%    G:                 inversion matrix
%    bQ:                data vector
%    name:              name of the figure
%    
% Output parameters:
%    image:             L-curves and cost functions
%    minimizeVector:    best solution minimizing both costs
%    infoVector:        informations related to minimization
%    minimizeValue:     regularization parameter

options             = IRset('RegParam', 'wgcv', 'NoStop', 'on',...
        'IterBar','off','verbosity','off');
[~, infoVector]     = IRhybrid_lsqr(G, bQ, options);
Reg                 = infoVector.StopReg;
It                  = Reg.It;
minimizeVector      = Reg.X;
minimizeValue       = Reg.RegP;
Res                 = infoVector.Rnrm;

image               =   figure('Name',name,...
    'NumberTitle','off','Position',[20,400,1200,1000],'visible','off');

plot(Res,'k','LineWidth',2);
hold on
scatter(It,Res(It),100,'filled','r');
text(It+10,Res(It),['Reg. Parameter = ' num2str(minimizeValue)],...
    'FontSize',15);
hold off
title('Regularization')
xlabel('Iteration');
ylabel('Cost Function');
SetFDefaults()

if outputLCurve == 1
    minimizeValue   = lCurveQc_k;
    options         = IRset('RegParam', lCurveQc_k, 'NoStop', 'on',...
        'IterBar','off','verbosity','off');
    [minimizeVector, infoVector] = IRhybrid_lsqr(G, bQ, options);
end


end