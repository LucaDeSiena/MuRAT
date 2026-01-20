function QcFrequency    =   Murat_imageQcFrequency(cf,...
    averageQcFrequency,sizeTitle,Qcf_title)
% function QcFrequency        =   Murat_imageQcFrequency(cf,...
%     averageQcFrequency,sizeTitle,Qcf_title)
%
% PLOTS the Qv vs Frequency relation and displays the fit
%
% Input parameters:
%    cf:                    central frequency
%    averageQcFrequency:    average inverse Qc
%    sizeTitle:             size of title font
%    Qcf_title:             title of peak delay figure
%
% Output parameters:
%    QcFrequency:           figure for Qc vs frequency

QcFrequency =   myfig(Qcf_title);

Qc_1        =   averageQcFrequency(1,:);
uncQc_1     =   averageQcFrequency(2,:);

d           =   log(Qc_1);
G           =   ones(length(cf),2);
G(:,2)      =   log(cf);

varQc       =   diag(uncQc_1.^(-2));

m           =   lsqlin(G,d);
covG        =   (G'*varQc*G)^(-1);
deltam      =   1.96*sqrt(diag(covG));

dpre        =   exp(G*m);
errorbar(cf,Qc_1,uncQc_1,'ko','LineWidth',2,'MarkerSize',12)

sQ0_1       =   exp(m(1));
unQ0_1      =   exp(deltam(1))*deltam(1);
unf0_1      =   deltam(2);

hold on
plot(cf,dpre,'k','LineWidth',2);
hold off
title(['Qc^{-1} = (' num2str(sQ0_1) '+/- ' num2str(unQ0_1) ') * f^{('...
    num2str(m(2)) ' +/- ' num2str(unf0_1) ')}'],'FontSize',sizeTitle);

xlabel('Frequency (Hz)','FontSize',sizeTitle,'FontWeight','bold',...
    'Color','k')
ylabel('Qc^{-1}','FontSize',sizeTitle,'FontWeight','bold','Color','k')
set(gca,'XScale','log','YScale','log')

SetFDefaults()