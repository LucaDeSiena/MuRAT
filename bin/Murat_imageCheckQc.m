function Qc_analysis        =   Murat_imageCheckQc(Qm_k,RZZ_k,...
    residualQc_k,luntot_Qc,Ac,sizeTitle,Qc_title,QcM)
% function Qc_analysis        =   Murat_imageCheckQc(Qm_k,RZZ_k,...
%     residualQc_k,luntot_Qc,Ac,sizeTitle,Qc_title,QcM)
%
% PLOTS the Qc checks
%
% Input parameters:
%    Qm_k:          inverse Qc data
%    RZZ_k:         uncertainty on Qc
%    luntot_Qc_k:   total lengths of rays for Qc
%    Ac:            inverse coda attenuation matrix
%    sizeTitle:     size of title font
%    Qc_title:      title of Qc figure
%    QcM:           Linearized or Non Linear Qc measurement
%
% Output parameters:
%    Qc_analysis:   figure for Qc check

Qc_analysis                 =   figure('Name',Qc_title,...
    'NumberTitle','off','Position',[20,400,1200,1000],'visible','off');
    
mQm                         =   mean(Qm_k);
Wc                          =   Murat_weighting(RZZ_k,QcM);
Gc                          =   Wc*Ac;
dck                         =   Wc*Qm_k;
[Uc,Sc,~]                   =   svd(Gc);
%%
% It starts with Qc relative to the ray length in the first plot.
subplot(2,1,1)
plot(luntot_Qc,Qm_k,'o','MarkerSize',6,'MarkerEdgeColor',[0 0 0])
hold on
plot(luntot_Qc,mQm*ones(length(luntot_Qc),1),'r-','LineWidth',2)
title('Dependence of Qc^{-1} on ray lengths');
xlabel('Ray length (km)','FontSize',sizeTitle,...
    'FontWeight','bold','Color','k')
ylabel('Qc^{-1}','FontSize',sizeTitle,...
    'FontWeight','bold','Color','k')
legend({'Qc^{-1}',cat(2,'<Qc> = ',num2str(1/mQm))},...
    'Location','northeast')
SetFDefaults()
%%
% A Picard plot evaluates the quality of the inversion.

subplot(2,1,2)
picard(Uc,diag(Sc),dck);
title(['Picard condition. Residual is: ' num2str(residualQc_k)]);
SetFDefaults()