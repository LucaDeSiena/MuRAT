function [image,parCond,para_map]    = ...
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
%    sizeTitle:         size of the title font
%
% Output parameters:
%    image:             image produced
%    para_condition:    condition on the parameters to image
%    para_map:          parameter map

para_map        =   Murat_unfoldXYZ(x,y,z);
condition       =   abs(modv_pd_k(:,4))>10^(-10);
parCond  =   para_map(condition,:);
pdCond    =   modv_pd_k(condition,4);
QcCond    =   modv_Qc_k(condition,4);

pdd             =   fitdist(pdCond,'ExtremeValue');

pdCond          =   pdCond - pdd.mu;
trepd           =   0.15*pdd.sigma;
mipdm           =   min(pdCond);
mapdm           =   max(pdCond);

Qcd             =   fitdist(QcCond,'GeneralizedExtremeValue');
QcCond    =   QcCond - Qcd.mu;
treQc           =   0.15*Qcd.sigma;
miQcm           =   min(QcCond);
maQcm           =   max(QcCond);

image           =...
    myfig('Parameter space separation');

% Precompute masks (each logical vector computed once)
m1      =   QcCond < -treQc & pdCond < -trepd;
m2      =   QcCond < -treQc & pdCond >  trepd;
m3      =   QcCond >  treQc & pdCond < -trepd;
m4      =   QcCond >  treQc & pdCond >  trepd;
m5      =  (QcCond>-treQc & QcCond<treQc) | (pdCond>-trepd & pdCond<trepd);

masks   =   {m1, m2, m3, m4, m5};
vals    =   [1, 2, 3, 4, 0];
sizes   =   [65, 65, 65, 65, 85];
colors  =   [ 0   0.8 0;   % green
            0   0.6 1;   % cyan
            1   0.6 0;   % orange
            1   0   0;   % red
            0.7 0.7 0.7];% gray

% Apply parCond assignments and scatter plots in a short loop
hold on
for k = 1:numel(masks)
    c   =   masks{k};
    if ~any(c), continue; end         % skip empty groups
    parCond(c,4)    =   vals(k);
    if k < 5
        scatter(QcCond(c), pdCond(c), sizes(k), 'filled', ...
                'MarkerFaceColor', colors(k,:));
    else
        scatter(QcCond(c), pdCond(c), sizes(k), 'filled', ...
                'MarkerFaceColor', colors(k,:), ...
                'MarkerEdgeColor', [1 1 1], 'LineWidth', 2);
    end
end

% Draw lines once (unchanged)
line([0 0], [mipdm mapdm], 'Color', [0 0 0], 'LineWidth', 3);
line([miQcm maQcm], [0 0], 'Color', [0 0 0], 'LineWidth', 3);

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
