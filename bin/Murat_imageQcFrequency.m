function QcFrequency        =   Murat_imageQcFrequency(cf,...
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

QcFrequency                 =   figure('Name',Qcf_title,...
    'NumberTitle','off','Position',[20,400,1200,1000],'visible','off');

Qc_1                        =   averageQcFrequency(1,:);
uncQc_1                     =   averageQcFrequency(2,:);

Qc                          =   averageQcFrequency(1,:).^(-1);
perQc                       =   uncQc_1./Qc_1;
uncQc                       =   perQc.*Qc;

d                           =   log(Qc);
G                           =   ones(length(cf),2);
G(:,2)                      =   log(cf);

varQc                       =   mean(perQc)^2;

m                           =   lsqlin(G,d);
covG                        =   varQc*(G'*G)^(-1);
deltam                      =   1.96*sqrt(diag(covG));

dpre                        =   exp(G*m);
errorbar(cf,Qc,uncQc,'ko','LineWidth',2,'MarkerSize',12)
hold on
plot(cf,dpre,'k','LineWidth',2);
hold off
title('Dependence of Qc on frequency','FontSize',sizeTitle);

xlabel('Frequency (Hz)','FontSize',sizeTitle,'FontWeight','bold',...
    'Color','k')
ylabel('Quality factor','FontSize',sizeTitle,...
    'FontWeight','bold','Color','k')
disp(['Constant = ' num2str(m(1)) ' +/- ' num2str(deltam(1))]);           
disp(['Frequency depencence = ', num2str(m(2)),...
    ' +/- ', num2str(deltam(2))]);

SetFDefaults()