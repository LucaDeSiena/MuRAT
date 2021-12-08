function [ParaMap,para_map]     = ...
    Murat_imageParametersMaps(par,para_map,x,y,z,X,Y,Z,evestaz_Qc,...
    sections,sizeTitle,FName_Parameters)
% function ParaMap                = ...
%     Murat_imageParametersMaps(par,para_map,x,y,z,X,Y,Z,evestaz_Qc,...
%     sections,sizeTitle,FName_Parameters)
%
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

for n = 1:length(par(:,1))
    locate_para                 =   para_map(:,1) == par(n,1) &...
        para_map(:,2) == par(n,2) & para_map(:,3) == par(n,3);
    para_map(locate_para,4)     =   par(n,4);
end
[~,~,~,para_3D]                 =   Murat_fold(x,y,z,para_map(:,4)-5);
un_X                            =   unique(para_3D);
fu                              =   find(un_X == -1, 1);
%%
% The four options are: (1) high for both (red); (2) low for both (green);
% (3) high for peak delays only (cyan); (4) high for inverse Qc only
% (orange).
if isempty(fu)
    
    cma_para                    =...
        [0.7 0.7 0.7;  0 0.8 0; 0 0.6 1; 1 0.6 0];
    HTick                       =...
        {'Average','Low Scattering\newline{Low Absorption}',...
        'High Scattering','High Absorption'};
    
else
    
    cma_para                    =...
        [0.7 0.7 0.7;  0 0.8 0; 0 0.6 1; 1 0.6 0; 1 0 0];
    HTick                       =...
        {'Average','Low Scattering\newline{Low Absorption}',...
        'High Scattering','High Absorption',...
        'High Scattering\newline{High Absorption}'};
    
end


ParaMap                         =...
    Murat_image3Dparameters(X,Y,Z,para_3D,...
    'bone',sections,evestaz_Qc,x,y,z,FName_Parameters);

colormap(cma_para)
colorbar('Ticks',un_X,'TickLabels',HTick);
title('Parameter separation map',...
    'FontSize',sizeTitle,'FontWeight','bold','Color','k');

end
