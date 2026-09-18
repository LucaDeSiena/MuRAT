function [lunpar, blocch, lunto, s, rayCrossing] =...
    Murat_segments(modv,rma)
% function [lunpar, blocch, lunto, s, rayCrossing] =...
%     Murat_segments(modv,rma)
%
% CREATES variables to compute ray-dependent inversion matrices and more
%
% Input parameters:
%    modv:          velocity model
%    rma:           ray in MuRAT format
%
% Output parameters:
%    lunpar:        the length of the ray in each crossed block
%    blocch:        crossed block numbers
%    lunto:         total ray length inside the model
%    s:             slowness in each crossed block
%    rayCrossing:   how many subsegments cross each block
%
% Structure:
% The ray traced by Murat_tracing is a polyline in columns 2:4 of rma.
% This routine splits each polyline segment at every crossed x/y/z cell
% boundary and accumulates the exact subsegment length in the containing
% cell. No source-receiver chord resampling is used.

nBlocks                     =   length(modv(:,1));
lunByBlock                  =   zeros(nBlocks,1);
rayCrossing                 =   zeros(nBlocks,1);

xv                          =   modv(:,1);
yv                          =   modv(:,2);
zv                          =   modv(:,3);
bv                          =   1:nBlocks;

lgx                         =   max(xv);
lgy                         =   max(yv);
lvg                         =   min(zv);
lgx1                        =   min(xv);
lgy1                        =   min(yv);
lvg1                        =   max(zv);

passox                      =   abs(find(xv~=xv(1),1,'first')-1);
passoy                      =   abs(find(yv~=yv(1),1,'first')-1);
passoz                      =   1;

deltastepx                  =   modv(1+passox,1)-modv(1,1);
deltastepy                  =   modv(1+passoy,2)-modv(1,2);
deltastepz                  =   modv(1+passoz,3)-modv(1,3);

mn                          =   [lgx-lgx1 lgy-lgy1 lvg-lvg1];
mn1                         =   max(mn);

BLx                         =   lgx1:deltastepx:lgx1+mn1;
BLy                         =   lgy1:deltastepy:lgy1+mn1;
BLv                         =   lvg1:deltastepz:lvg;

fallr                       =   find(rma(:,2)>0 & rma(:,3)>0);
if length(fallr) < 2
    blocch                  =   zeros(0,1);
    lunpar                  =   zeros(0,1);
    lunto                   =   0;
    s                       =   zeros(0,1);
    return
end

rayPoints                   =   rma(fallr,2:4);
rayStep                     =   sqrt(sum(diff(rayPoints,1,1).^2,2));
keepRayPoint                =   [true; rayStep>0];
rayPoints                   =   rayPoints(keepRayPoint,:);

for ir = 1:size(rayPoints,1)-1
    p0                      =   rayPoints(ir,:);
    p1                      =   rayPoints(ir+1,:);
    dp                      =   p1-p0;
    segLength               =   sqrt(sum(dp.^2));
    if segLength == 0
        continue
    end

    tBreak                  =   [0 1];

    if dp(1) ~= 0
        tx                  =   (BLx-p0(1))/dp(1);
        tBreak              =   [tBreak tx(tx>0 & tx<1)];
    end
    if dp(2) ~= 0
        ty                  =   (BLy-p0(2))/dp(2);
        tBreak              =   [tBreak ty(ty>0 & ty<1)];
    end
    if dp(3) ~= 0
        tz                  =   (BLv-p0(3))/dp(3);
        tBreak              =   [tBreak tz(tz>0 & tz<1)];
    end

    tBreak                  =   sort(tBreak);
    tBreak                  =   tBreak([true diff(tBreak)>1e-10]);

    for it = 1:length(tBreak)-1
        t0                  =   tBreak(it);
        t1                  =   tBreak(it+1);
        if t1 <= t0
            continue
        end

        tMid                =   (t0+t1)/2;
        pMid                =   p0+tMid*dp;
        subLength           =   segLength*(t1-t0);

        inBlock             =   xv <= pMid(1) & pMid(1) < xv+deltastepx...
            & yv <= pMid(2) & pMid(2) < yv+deltastepy...
            & zv+deltastepz <= pMid(3) & pMid(3) < zv;
        block               =   bv(inBlock>0)+1;

        if isempty(block)
            continue
        end

        block               =   block(1);
        if block > nBlocks || block < 1
            continue
        end

        lunByBlock(block)   =   lunByBlock(block)+subLength;
        rayCrossing(block)  =   rayCrossing(block)+1;
    end
end

blocch                      =   find(lunByBlock>0);
lunpar                      =   lunByBlock(blocch);
lunto                       =   sum(lunpar);
s                           =   zeros(size(lunpar));

if ~isempty(blocch)
    s                       =   1./modv(blocch,4);
end
end
