function [image,para_condition,para_map]= ...
    Murat_imageParameters(x,y,z,modv_pd_k,modv_Qc_k,sizeTitle)
% function [image,para_condition,para_map]= ...
%     Murat_imageParameters(x,y,z,modv_pd_k,modv_Qc_k,sizeTitle)
% PLOTS the parameter models in their space and on a 3D map
%
% Input parameters:
%    x:                 x vector
%    y:                 y vector
%    z:                 z vector
%    modv_pd_k:         results of peak delay imaging
%    modv_Qc_k:         results of Qc imaging
%
% Output parameters:
%    image:             image produced
%    para_condition:    condition on the parameters to image
%    para_map:          parameter map

para_map                                =   Murat_unfoldXYZ(x,y,z);
condition                               =   abs(modv_pd_k(:,4))>10^(-10);
para_condition                          =   para_map(condition,:);
pd_condition                            =   modv_pd_k(condition,4);
Qc_condition                            =   modv_Qc_k(condition,4);

pdd                                     =...
    fitdist(pd_condition,'ExtremeValue');

pd_condition                            =   pd_condition - pdd.mu;
trepd                                   =   0.15*pdd.sigma;
mipdm                                   =   min(pd_condition);
mapdm                                   =   max(pd_condition);

Qcd                                     =...
    fitdist(Qc_condition,'GeneralizedExtremeValue');
Qc_condition                            =   Qc_condition - Qcd.mu;
treQc                                   =   0.15*Qcd.sigma;
miQcm                                   =   min(Qc_condition);
maQcm                                   =   max(Qc_condition);

image                                   =...
    figure('Name','Parameter space separation',...
    'NumberTitle','off','Position',[300,200,1200,1000],'visible','off');

c                                       =...
    Qc_condition<-treQc & pd_condition<-trepd;
para_condition(c,4)                     =   1;
scatter(Qc_condition(c),pd_condition(c),65,'filled',...
    'MarkerFaceColor',[0 0.8 0])

hold on
line([0 0],[mipdm mapdm],'Color',[0 0 0],'LineWidth',3)
line([miQcm maQcm],[0 0],'Color',[0 0 0],'LineWidth',3)

c                                       =...
    Qc_condition<-treQc & pd_condition>trepd;
para_condition(c,4)                     =   2;
scatter(Qc_condition(c),pd_condition(c),65,'filled',...
    'MarkerFaceColor',[0 0.6 1])

c                                       =...
    Qc_condition>treQc & pd_condition<-trepd;
para_condition(c,4)                     =   3;
scatter(Qc_condition(c),pd_condition(c),65,'filled',...
    'MarkerFaceColor',[1 0.6 0])

c                                       =...
    Qc_condition>treQc & pd_condition>trepd;
para_condition(c,4)                     =   4;
scatter(Qc_condition(c),pd_condition(c),65,'filled',...
    'MarkerFaceColor',[1 0 0])

c                                       =...
    (Qc_condition>-treQc & Qc_condition<treQc) |...
    (pd_condition>-trepd & pd_condition<trepd);
para_condition(c,4)                     =   0;
scatter(Qc_condition(c),pd_condition(c),85,'filled','MarkerFaceColor',...
    [0.7 0.7 0.7],'MarkerEdgeColor',[1 1 1],'LineWidth',2)

xlim([miQcm maQcm]);
ylim([mipdm mapdm]);

SetFDefaults();

xlabel('Qc','FontSize',sizeTitle,'FontWeight','bold','Color','k')
ylabel('Log. peak delay','FontSize',sizeTitle,...
    'FontWeight','bold','Color','k')
title('Parameter space plot',...
    'FontSize',sizeTitle,'FontWeight','bold','Color','k');
axis square

hold off
