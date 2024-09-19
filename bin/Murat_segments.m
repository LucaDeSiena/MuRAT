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
%    lunpar:        the lenght of the segments in each block
%    blocch:        the corresponding block number
%    lunto:         the total length of the segments
%    s:             the slowness in each block of the vel. model
%    rayCrossing:   how many rays cross each block
%             
% Structure:
% The velocity model is a x,y,z grid with regular step IN METERS, and depth
%   is positive. This means that the third coordinate of the input velocity
%   model must be flipped, and altitude become negative.
%   The variable "rma" is the ray traced in the velocity model by the
%   "tracing" routine. This variable has already the right orientation.
% The code takes two points of the ray and checks if a
%   boundary is crossed. If not, the length of the segment is added to
%   the previous part. It the bounday is crossed, the segment is broken,
%   and a new segment starts. The ray gets interpolated into smaller
%   segments for hit check.

lunpar                      =	zeros(100,1);
blocch                      =	zeros(100,1);
rayCrossing                 =	zeros(length(modv(:,1)),1);

lgx                         =   max(modv(:,1));
lgy                         =   max(modv(:,2));
lvg                         =   min(modv(:,3));
lgx1                        =   min(modv(:,1));
lgy1                        =   min(modv(:,2));
lvg1                        =   max(modv(:,3));

passox                      =   abs(find(modv(:,1)~=modv(1,1),1,'first')-1);
passoy                      =   abs(find(modv(:,2)~=modv(1,2),1,'first')-1);
passoz                      =   1;

% Find the index of variations in the velocity model
deltastepx                  =   modv(1+passox,1)-modv(1,1);
deltastepy                  =   modv(1+passoy,2)-modv(1,2);

% Decreases altitude
deltastepz                  =   modv(1+passoz,3)-modv(1,3);

% Creation of the grid, having constant steps
mn                          =   [lgx-lgx1 lgy-lgy1 lvg-lvg1];

% Find the longest between x, y ,z
mn1                         =   max(mn);

% Build up the corresponding vectors
BLx                         =   lgx1:deltastepx:lgx1+mn1;
BLy                         =   lgy1:deltastepy:lgy1+mn1;
BLv                         =   lvg1:deltastepz:lvg;

% Check that all x-y elements are more than zero
fallr                       =   find(rma(:,2)>0 & rma(:,3)>0);

% Define the interpolated vectors for hit count
xpolo_noint                 =   rma(fallr,2);
ypolo_noint                 =   rma(fallr,3);
zpolo_noint                 =   rma(fallr,4);

tot                         =   1000;
xpolo                       =...
    linspace(xpolo_noint(1),xpolo_noint(end),tot);
ypolo                       =...
    linspace(ypolo_noint(1),ypolo_noint(end),tot);
zpolo                       =...
    linspace(zpolo_noint(1),zpolo_noint(end),tot);

% Read hypocenter coordinates
sorg                        =   [xpolo(1) ypolo(1) zpolo(1)];

% Read station coordinates
ricev                       =   [xpolo(end) ypolo(end) zpolo(end)];

index                       =   1;
xv                          =   modv(:,1);
yv                          =   modv(:,2);
zv                          =   modv(:,3);

% Number corresponding to each block of the grid
bv                          =   1:length(xv);
    
%%
% Find the block corresponding to the source
bS                          =   xv <= sorg(1) & sorg(1) < xv+deltastepx...
    & yv <= sorg(2) & sorg(2) < yv+deltastepy...
    & zv+deltastepz <= sorg(3) & sorg(3) < zv;

% This is the block, added +1 as reference is the deepest point to the SW
bSS                         =   bv(bS>0)+1;

if isempty(bSS)
    error('Source outside velocity model!!')
end
%%    
% Creation of file inte, which samples the crossing points of the ray 
inte                        =   zeros(100,6);
inte(1,1:6)                 =   [sorg(1),sorg(2),sorg(3),0,0,0];
inte(2,6)                   =   bSS;

diffLoc                     =   zeros(tot,1);
k1                          =   1;
for k = 2:tot
    ewS                     =	xpolo(k-1); 
    nsS                     =   ypolo(k-1);
    vS                      =   zpolo(k-1);
    ewR                     =   xpolo(k);
    nsR                     =   ypolo(k);
    vR                      =   zpolo(k);

    diffLoc(k)              =   sqrt((ewR-ewS)^2+(nsR-nsS)^2+(vR-vS)^2);
    
    cx                      =   (ewS <= BLx & ewR > BLx) |...
        (ewS >=BLx & ewR < BLx);
    [a_x, b_x]              = find(cx);
    cy                      =   (nsS <= BLy & nsR > BLy) |...
        (nsS >=BLy & nsR < BLy);
    [a_y, b_y]              = find(cy);
    cv                      =   (vS <= BLv & vR > BLv) |...
        (vS >= BLv & vR < BLv);
    [a_z, b_z]              = find(cv);
    
    % x crossing
    if  a_x > 0
        index               =   index + 1;
        inte(index,1:5)     =   [BLx(b_x),nsR,vR,index-1,0];
        lung                =   sum(diffLoc(k1:k));
        k1                  =   k;
        inte(index,5)       =   lung;
        if ewS <= ewR
            if inte(index,6)+passox < max(bv)
                inte(index+1,6) = inte(index,6)+passox;
            else
                break
            end
        elseif ewS > ewR
            if inte(index,6)-passox > 0
                inte(index+1,6) = inte(index,6)-passox;
            else
                break
            end
        end
    end

    % y crossing
    if a_y > 0
        index               =   index + 1;
        inte(index,1:5)     =   [ewR,BLy(b_y),vR,index-1,0];
        lung                =   sum(diffLoc(k1:k));
        k1                  =   k;
        inte(index,5)       =   lung;
        if nsS <= nsR
            if inte(index,6)+passox < max(bv)
                inte(index+1,6) =   inte(index,6)+passoy;
            else
                break
            end
        elseif nsS > nsR
            if inte(index,6)-passoy > 0
                inte(index+1,6) =   inte(index,6)-passoy;
            else
                break
            end
        end
    end

    % z crossing
    if a_z > 0
        index               =   index + 1;
        inte(index,1:5)     =   [ewR,nsR,BLv(b_z),index-1,0];
        lung                =   sum(diffLoc(k1:k));
        k1                  =   k;
        inte(index,5)       =   lung;
        if vS <= vR
            if inte(index,6)+passoz < max(bv)
                inte(index+1,6) = inte(index,6)-passoz;
            else
                break
            end
        elseif vS > vR
            if inte(index,6)+passoz > 0
                inte(index+1,6) = inte(index,6)+passoz;
            else
                break
            end
        end
        
    end
    continue
end

% In case the receiver is in the grid.
if k == tot

    % receiver location
    bR                  =   xv < ricev(1) & ricev(1)<xv+deltastepx...
        & yv < ricev(2) & ricev(2)<yv+deltastepy...
        & zv+deltastepz<ricev(3) & ricev(3)<zv;
    bRR                 =   bv(bR>0);

    if isempty(bRR) == 0
        index           =   index + 1;
        lung            =   sum(diffLoc(k1:end));
        inte(index,1:5) =   [ricev(1),ricev(2),ricev(3),index-1,lung];
    end
end

% Variable inte monitors the points where the rays cross the grid
linte                   =   length(inte(:,1));

lunpar(1:linte,1)       =   inte(:,5);
blocch(1:linte,1)       =   inte(:,6);
lunto(1,1)              =   sum(lunpar(:,1));

s                       =   zeros(100,1);
no                      =   blocch>length(modv(:,1));
blocch(no)              =   0;
lunpar(no)              =   0;
bff                     =   find(blocch);

if bff > 0
    s(bff)              =   1./modv(blocch(bff),4);
    rayCrossing(blocch(bff))=   rayCrossing(blocch(bff))+1;
end
    
end
