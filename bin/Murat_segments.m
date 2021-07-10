function [lunpar, blocch, lunto, s, rayCrossing] =...
    Murat_segments(modv,rma)
%   segments  Creates files with
%             1. the total length of the segments (lunto)
%             2. the slowness (s) in each block of the vel. model
%             3. the lenght of the segments (lunpar) in each block
%             4. the corresponding file containing the number of blocks
%             crossed
%             
%             	
%==========================================================================
% The velocity model is a x,y,z grid with regular step IN METERS, and depth
% is positive. This means that the third coordinate of the input velocity
% model must be flipped, and altitude become negative.
% The variable "rma" is the ray traced in the velocity model by the
% "tracing" routine. This variable has already the right orientation.
%==========================================================================
% Set up the variables to save
lunpar                  = zeros(100,1); % Partial lengths in the blocks
blocch                  = zeros(100,1); % Numbers corresponding to the blocks
rayCrossing             = zeros(length(modv(:,1)),1); %Crosses ray

lgx                     = max(modv(:,1));%x maximum
lgy                     = max(modv(:,2));%y maximum
lvg                     = min(modv(:,3));% maximum depth
lgx1                    = min(modv(:,1));%x minimum
lgy1                    = min(modv(:,2));%y minimum
lvg1                    = max(modv(:,3));% minimum depth

% Check steps in the three directions, step for the z-coordinate is 1
passox                  = find(modv(:,1)~=modv(1,1),1,'first')-1;
passoy                  = find(modv(:,2)~=modv(1,2),1,'first')-1;
passoz                  = 1;

% Find the index of variations in the velocity model
deltastepx              = modv(1+passox,1);
deltastepy              = modv(1+passoy,2);

% Crucial to decrease altitude
deltastepz              = modv(1+passoz,3)-modv(1,3);

% Creation of the grid, having constant steps
mn                      = [lgx-lgx1 lgy-lgy1 lvg-lvg1];

% Find the longest  between x, y ,z
mn1                     = max(mn);

% Build up the corresponding vectors
BLx                     = lgx1:deltastepx:lgx1+mn1;
BLy                     = lgy1:deltastepy:lgy1+mn1;
BLv                     = lvg1:deltastepz:lvg;

% Check that all x-y elements are more than zero
fallr                   = find(rma(:,1)>0 & rma(:,2)>0);

% Define the vectors
xpolo                   = rma(fallr,2);
ypolo                   = rma(fallr,3);
% Need to switch as rma is with positive depths
zpolo                   = rma(fallr,4);
tot                     = length(xpolo);

sorg = [xpolo(1) ypolo(1) zpolo(1)];% read hypocenter coordinates
ricev = [xpolo(end) ypolo(end) zpolo(end)];% read station coordinates

index = 1;
xv = modv(:,1); %WE coordinate
yv = modv(:,2); %SN coordinate
zv = modv(:,3); %Z coordinate
bv = 1:length(xv);% number corresponding to each block of the grid
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the block corresponding to the source
bS = xv<=sorg(1) & sorg(1)<xv+deltastepx & yv<=sorg(2) &...
    sorg(2)<yv+deltastepy  & zv+deltastepz<=sorg(3) & sorg(3)<zv;

% This is the block, added +1 as the reference is the deepest point to the
% South-West
bSS = bv(bS>0)+1;

if isempty(bSS)
    error('Source outside velocity model!!')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%creation of file inte, which samples the crossing points of the ray 
inte = zeros(100,6);
inte(1,1:6) = [sorg(1),sorg(2),sorg(3),0,0,0]; %start at the source
inte(2,6) = bSS;

% Loop to find the length of the segments of the ray in each block of
% the grid. The code takes two points of the ray and checks if a
% boundary is crossed. If not, the length of the segment is added to
% the previous part. It the bounday is crossed, the segment is broken,
% and a new segment starts.
dif = zeros(tot,1);
k1= 1;
for k = 2:tot
    ewS = xpolo(k-1); %first WE coord. 
    nsS = ypolo(k-1); %first SN coord.
    vS =  zpolo(k-1); %first vertical coord.
    ewR = xpolo(k); %second WE coord.
    nsR = ypolo(k); %second SN coord.
    vR =  zpolo(k); %second vertical coord.

    dif(k)=sqrt((ewR-ewS)^2+(nsR-nsS)^2+(vR-vS)^2); % length of each segment
    
    cx = (ewS <= BLx & ewR > BLx) | (ewS >=BLx & ewR < BLx); %crossing WE
    cy = (nsS <= BLy & nsR > BLy) | (nsS >=BLy & nsR < BLy); % crossing SN
    cv = (vS <= BLv & vR > BLv) | (vS >= BLv & vR < BLv); %crosses vertical
    % x crossing
    if  find(cx) > 0
        n = find(cx>0, 1, 'first');
        index = index + 1;
        inte(index,1:5) = [BLx(n),nsR,vR,index-1,0];
        lung = sum(dif(k1:k));
        k1=k;
        inte(index,5) = lung;
        if ewS <= ewR
            if inte(index,6)+passox < max(bv)
                inte(index+1,6) = inte(index,6)+passox;
            else
                break %in case the grid is left by the ray
            end
        elseif ewS > ewR
            if inte(index,6)-passox > 0
                inte(index+1,6) = inte(index,6)-passox;
            else
                break %in case the grid is left by the ray
            end
        end
    end      
        
    % y crossing    
    if find(cy) > 0
        n = find(cy>0, 1, 'first');
        index = index + 1;
        inte(index,1:5) = [ewR,BLy(n),vR,index-1,0];
        lung = sum(dif(k1:k));
        k1=k;
        inte(index,5) = lung;
        if nsS <= nsR
            if inte(index,6)+passox < max(bv)
                inte(index+1,6) = inte(index,6)+passoy;
            else
                break %in case the grid is left by the ray
            end
        elseif nsS > nsR
            if inte(index,6)-passoy > 0
                inte(index+1,6) = inte(index,6)-passoy;
            else
                break %in case the grid is left by the ray
            end
        end
    end   
        
    % z crossing    
    if find(cv) > 0
        n = find(cv>0, 1, 'first');
        index = index + 1;
        inte(index,1:5) = [ewR,nsR,BLv(n),index-1,0];
        lung = sum(dif(k1:k));
        inte(index,5) = lung;
        if vS <= vR
            if inte(index,6)+passoz < max(bv)
                inte(index+1,6) = inte(index,6)-passoz;
            else
                break %in case the grid is left by the ray
            end
        elseif vS > vR
            if inte(index,6)+passoz > 0
                inte(index+1,6) = inte(index,6)+passoz;
            else
                break %in case the grid is left by the ray
            end
        end
       
    end
    continue
end
    
% In case the receiver is in the grid.
if k == tot
    % receiver location
    bR = xv < ricev(1) & ricev(1)<xv+deltastepx...
        & yv < ricev(2) & ricev(2)<yv+deltastepy...
        & zv+deltastepz<ricev(3) & ricev(3)<zv;
    bRR = bv(bR>0);

    if isempty(bRR)==0
        index = index + 1;
        lung = sum(dif(k1:end));
        inte(index,1:5) = [ricev(1),ricev(2),ricev(3),index-1,lung];
    end
end

% Variable inte monitors the points where the rays cross the grid
linte = length(inte(:,1));
% Variable lunpar contains the lengths of each segment in each block
lunpar(1:linte,1) = inte(:,5);
% Variable blocch relates the length of the segment in lunpar to the
% corresponding block of the grid
blocch(1:linte,1) = inte(:,6);
% Total length of each ray
lunto(1,1) = sum(lunpar(:,1));


%==========================================================================
% The following loop creates a variable containing the slowness from the
% velocity model for the blocks crossed by each ray (s) and adds to the
% velocity model a sixth column containing the number of rays crossing each
% block.
%==========================================================================

s = zeros(100,1);
no=blocch>length(modv(:,1));
blocch(no)=0;
lunpar(no)=0;
bff=find(blocch);

if bff > 0
    s(bff)=1./modv(blocch(bff),4);
    rayCrossing(blocch(bff)) = rayCrossing(blocch(bff))+1;
end

end
