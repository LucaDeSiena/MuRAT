function [Apd_i, A_i, luntot_i, rma_i, rayCrossing_i]  =...
    Murat_rays(modv,gridD,pvel,sst_i)
% function [Apd_i, A_i, luntot_i, rma_i, rayCrossing_i]  =...
%     Murat_rays(modv,gridD,pvel,sst_i)
%
% RAY TRACING and computation of SEGMENTS inside each cell.
%
% Input parameters:
%    modv:          velocity model for inversion
%    gridD:         grid for ray tracing
%    pvel:          velocity model for ray tracing
%    sst_i:         info on source and stations
%
% Output parameters:
%    Apd_i, Ac_i:   inversion matrices for peak delay and coda
%    luntot_i:      total length of ray
%    rma_i:         rays in the three components
%    rayCrossing_i: hit counts

nCells  =   size(modv,1);

Apd_i   =   zeros(1,nCells);
A_i     =   zeros(1,nCells);

ray     =   reshape(sst_i,3,2);
rma     =   Murat_tracing(ray,gridD,pvel);

%Calculation of segment lengths
[lunpar, blocch,luntot_i, s, rayCrossing_i] =   Murat_segments(modv,rma);

hitMask         =   blocch > 0;
cellsHit        =   blocch(hitMask);

Apd_i(cellsHit) =   lunpar(hitMask)/1000;
A_i(cellsHit)   =   -lunpar(hitMask)/1000.*s(hitMask);

% Interpolation for plotting rays
lrma            =   rma(end,1);
idx             =   round(linspace(1, max(1, lrma), 100));
rma_i           =   rma(idx, :);
end
