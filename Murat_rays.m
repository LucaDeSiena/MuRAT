function [Apd_i, A_i, luntot_i, rma_i,rayCrossing_i]  =...
    Murat_rays(Murat,sst_i)

modv                                    =   Murat.input.modv;
gridD                                   =   Murat.input.gridD;
pvel                                    =   Murat.input.pvel;

Apd_i                                   =   zeros(1,length(modv(:,1)));
A_i                                     =   zeros(1,length(modv(:,1)));

%RAY TRACING
ray                                     =   reshape(sst_i,3,2);
rma                                     =   Murat_tracing(ray,gridD,pvel);

%CALCULATON OF SEGMENT LENGTHS
[lunpar, blocch,...
    luntot_i, s, rayCrossing_i]         =   Murat_segments(modv,rma);

lb                                      =   blocch>0;
Apd_i(blocch(lb))                       =   lunpar(lb)/1000;
A_i(blocch(lb))                         =   -lunpar(lb)/1000.*s(lb);

% Interpolation for plotting rays
lrma                                    =   rma(end,1);
in_rma                                  =   floor(linspace(1,lrma,100));
rma_i                                   =   rma(in_rma,:);
