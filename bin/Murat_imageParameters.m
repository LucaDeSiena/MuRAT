function [image,par,para_map]    = ...
    Murat_imageParameters(x,y,z,modv1,modv2,visib)
%PLOTS the parameter models in their space and on a 3D map

para_map                         =   Murat_unfold(x',y',z');
condition                        =   abs(modv1(:,4))>10^(-10);
xyz_condition                    =   para_map(condition,:);
pd_condition                     =   modv1(condition,4);
pdef                             =   [xyz_condition pd_condition];
Qc_condition                     =   modv2(condition,4);
Qcef                             =   [xyz_condition Qc_condition];

pdps                             =   pdef(:,4);
Qps                              =   Qcef(:,4);

pdd                              =   fitdist(pdps,'ExtremeValue');
% figure
% mipdm                   =   min(pdef(:,4));
% mapdm                   =   max(pdef(:,4));
% miQcm                   =   min(Qcef(:,4));
% maQcm                   =   max(Qcef(:,4));
% P = histogram(pdps,'Normalization','pdf','FaceColor',[.9 .9 .9]);
% ylabel('Probability Density');
% xlabel('Peak Delay');
% 
% xgrid = linspace(mipdm,mapdm,1000)';
% pdfEst = pdf(pdd,xgrid);
% line(xgrid,pdfEst)

pdps                            =   pdps - pdd.mu;
trepd                           =   0.01*pdd.sigma;
mipdm                           =   min(pdps);
mapdm                           =   max(pdps);

Qcd                             =   fitdist(Qps,'GeneralizedExtremeValue');
% figure
% Q = histogram(Qps,'Normalization','pdf','FaceColor',[.9 .9 .9]);
% ylabel('Probability Density');
% xlabel('Qc');
% 
% xgrid = linspace(miQcm,maQcm,1000)';
% pdfEst = pdf(Qcd,xgrid);
% line(xgrid,pdfEst)

Qps                             =   Qps - Qcd.mu;
treQc                           =   0.01*Qcd.sigma;
miQcm                           =   min(Qps);
maQcm                           =   max(Qps);

image                           =...
    figure('Name','Parameter space separation',...
    'NumberTitle','off','visible',visib,'Position',[300,200,1200,1000]);

par                             =   pdef(:,1:3);

c                               =   Qps<-treQc & pdps<-trepd;
par(c,4)                        =   1;
scatter(Qps(c),pdps(c),65,'filled','MarkerFaceColor',[0 0.8 0])

hold on
line([0 0],[mipdm mapdm],'Color',[0 0 0],...
    'LineWidth',3)
line([miQcm maQcm],[0 0],'Color',[0 0 0],...
    'LineWidth',3)
c                               =   Qps<-treQc & pdps>trepd;
par(c,4)                        =   2;
scatter(Qps(c),pdps(c),65,'filled','MarkerFaceColor',[0 0.6 1])

c                               =   Qps>treQc & pdps<-trepd;
par(c,4)                        =   3;
scatter(Qps(c),pdps(c),65,'filled','MarkerFaceColor',[1 0.6 0])

c                               =   Qps>treQc & pdps>trepd;
par(c,4)                        =   4;
scatter(Qps(c),pdps(c),65,'filled','MarkerFaceColor',[1 0 0])

c                               =   (Qps>-treQc & Qps<treQc) |...
    (pdps>-trepd & pdps<trepd);
par(c,4)                        =   0;
scatter(Qps(c),pdps(c),85,'filled','MarkerFaceColor',[0.7 0.7 0.7],...
    'MarkerEdgeColor',[1 1 1],'LineWidth',2)
xlim([miQcm maQcm]);
ylim([mipdm mapdm]);

SetFDefaults();
hold off
