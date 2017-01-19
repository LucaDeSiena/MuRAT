function [lunparz, blocchi, luntot, inblocchi, sb] = segments2(modv,originz,um)


%   segments  Creates files with
%             1. the total length of the segments (luntot)
%             2. the slowness (sb) in each block of the vel. model
%             3. the lenght of the segments (lunparz) in each block
%             4. the corresponding file containing the number of blocks
%             crossed
%             5. the coordinates of the blocks effectively crossed with the number of rays crossing them (inblocchi)
%
%             	
%==========================================================================
% The ray-path contained in the folder "rays" are used.
% The files are ascii files containing the cartesian coordinates
% of the rays in the reference system of the velocity model, starting at
% the source and with positive depth. The velocity model is a x,y,z grid
% with regular step IN METERS
%==========================================================================
lgx = max(modv(:,1));%x maximum
lgy = max(modv(:,2));%y maximum
lvg = -min(modv(:,3));% maximum depth
lgx1 = min(modv(:,1));%x minimum
lgy1 = min(modv(:,2));%y minimum
lvg1 = -max(modv(:,3));% minimum depth

% Check steps in the three directions, step for the z-coordinate is 1
passox = find(modv(:,1)~=modv(1,1),1,'first')-1;
passoy = find(modv(:,2)~=modv(1,2),1,'first')-1;
passoz = 1;
deltastepx = modv(1+passox,1)-modv(1,1);
deltastepy = modv(1+passoy,2)-modv(1,2);
deltastepz = modv(1+passoz,3)-modv(1,3);

%creation of the grid, having constant steps
mn = [lgx-lgx1 lgy-lgy1 lvg-lvg1];
mn1 = max(mn);
BLx = lgx1:deltastepx:lgx1+mn1;
BLy = lgy1:deltastepy:lgy1+mn1;
BLv = lvg:deltastepz:lvg1;

% List of the paths of the ray-files - in rays
list=dir('./rays');  %get info of files/folders in current directory
isfile=~[list.isdir]; %determine index of files vs folders
filenames={list(isfile).name}; %create cell array of file names
listaraggi = filenames';
for i = 1:length(listaraggi)
    listaraggi{i,1}=cat(2,'./rays/',listaraggi{i});
end

%==========================================================================
% A loop starts for the rays in the first velocity model, q is the index
% referred to each ray file. The order of the ray files must be the same of
% the traces and time picking files.
%
%EXAMPLE: 2 rays recorded at 1 three component station.
%         list file for the rays: 6 rows
%         list file for the traces: 6 rows (1 for each component)
%         time picking file: 6 rows (1 for each component even if
%         the pickings are identical)
%==========================================================================
lr = length(listaraggi(:,1));
luntot = zeros(1,lr);
lunparz = zeros(100,lr);
blocchi = zeros(100,lr);
indexr = 0;

for q = 1:lr; %for loop, the index corresponds to the ray
    indexr = indexr+1;
    display(indexr);
    fileID=fopen(listaraggi{q,1});
    allr=textscan(fileID,'%f %f %f %f %f');% read ray
    fclose(fileID);
    xpolo=[allr{2}];
    ypolo=[allr{3}];
    zpolo=[allr{4}];
    tot = length(xpolo);

    %the coordinate of the ray with topography - depending on ray-tracing
    %remember that the ray coordinate must be positive - opposite to the
    %velocity model. See the MSH example.
    zpolo = zpolo-originz;
    sorg = [xpolo(1) ypolo(1) -zpolo(1)];% read hypocenter coordinates
    ricev = [xpolo(end) ypolo(end) -zpolo(end)];% read station coordinates
    index = 1;
    xv = modv(:,1); %WE coordinate
    yv = modv(:,2); %SN coordinate
    zv = modv(:,3); %Z coordinate
    bv = 1:length(xv);% number corresponding to each block of the grid
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bS = xv<=sorg(1) & sorg(1)<xv+deltastepx & yv<=sorg(2) &...
        sorg(2)<yv+deltastepy  & zv<=sorg(3) & sorg(3)<zv-deltastepz;
    bSS = bv(bS>0);
    non = 1;
    while isempty(bSS) && non<tot
        non=non+1;
        sorg = [xpolo(non) ypolo(non) -zpolo(non)];% read hypocenter coordinates
        bS = xv<=sorg(1) & sorg(1)<(xv+deltastepx) & yv<=sorg(2) &...
            sorg(2)<yv+deltastepy  & zv<=sorg(3) & sorg(3)<zv-deltastepz;
        bSS = bv(bS>0);
    end
    if isempty(bSS)
        continue
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %creation of file inte, which samples the crossing points of the ray 
    inte = zeros(100,6);
    inte(1,1:6) = [sorg(1),sorg(2),sorg(3),0,0,0]; %start at the source
    inte(2,6) = bSS;
    
    % Loop to find the length of the segments of the ray in each block of
    % the grid. The code takes two ponts of the ray and checks if a
    % boundary is crossed. If not, the length of the segment is added to
    % the previous part. It the bounday is crossed, the segment is broken,
    % and a new segment starts.
    dif = zeros(tot,1);
    k1= 1;
    for k = non+1:tot
        ewS = xpolo(k-1); %first WE coord. 
        nsS = ypolo(k-1); %first SN coord.
        vS =  zpolo(k-1); %first vertical coord.
        ewR = xpolo(k); %second WE coord.
        nsR = ypolo(k); %second SN coord.
        vR =  zpolo(k); %second vertical coord.

        dif(k)=sqrt((ewR-ewS)^2+(nsR-nsS)^2+(vR-vS)^2); % length of each segment
        
        cx = (ewS <= BLx & ewR > BLx) | (ewS >=BLx & ewR < BLx); %crossing WE
        cy = (nsS <= BLy & nsR > BLy) | (nsS >=BLy & nsR < BLy); % crossing SN
        cv = (vS <= BLv & vR > BLv) | (vS >= BLv & vR < BLv); %crossing vertical
        % x crossing
        if  find(cx) > 0
            n = find(cx>0, 1, 'first');
            index = index + 1;
            inte(index,1:5) = [BLx(n),nsR,-vR,index-1,0];
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
            inte(index,1:5) = [ewR,BLy(n),-vR,index-1,0];
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
            n = cv>0;
            index = index + 1;
            inte(index,1:5) = [ewR,nsR,-BLv(n),index-1,0];
            lung = sum(dif(k1:k));
            inte(index,5) = lung;
            if vS <= vR
                if inte(index,6)+passoz < max(bv)
                    inte(index+1,6) = inte(index,6)+passoz;
                else
                    break %in case the grid is left by the ray
                end
            elseif vS > vR
                if inte(index,6)+passoz > 0
                    inte(index+1,6) = inte(index,6)-passoz;
                else
                    break %in case the grid is left by the ray
                end
            end
           
        end
        continue
    end
    
    % In case the receiver is in the grid.
    if non == tot
        % receiver location
        bR = xv < ricev(1) & ricev(1)<xv+deltastepx...
            & yv < ricev(2) & ricev(2)<yv+deltastepy...
            & zv<ricev(3) & ricev(3)<zv-deltastepz;
        bRR = bv(bR>0);
    
        if isempty(bRR)==0
            index = index + 1;
            lung = sum(dif(k1:end));
            inte(end,1:5) = [ricev(1),ricev(2),ricev(3),index-1,lung];
        end
    end
    % Variable inte monitors the points where the rays cross the grid
    linte = length(inte(:,1));
    % Variable lunparz contains the lengths of each segment in each block
    lunparz(1:linte,indexr) = inte(:,5);
    % Variable blocchi relates the length of the segment in lunparz to the
    % corresponding block of the grid
    blocchi(1:linte,indexr) = inte(:,6);
    % Total length of each ray
    luntot(1,indexr) = sum(lunparz(:,indexr));
    %Theoretical travel-time of the ray
end

luntot = luntot/um;
lunparz = lunparz/um;

%==========================================================================
% The following loop creates a variable containing the slowness from the
% velocity model for the blocks crossed by each ray (sb) and adds to the
% velocity model a sixth column containing the number of rays crossing each
% block.
%==========================================================================

sb = zeros(100,lr);
count = (zeros(1,length(modv(:,1))));
no=blocchi>length(modv(:,1));
blocchi(no)=0;
lunparz(no)=0;
for j = 1:lr
    bff=find(blocchi(:,j));
    if bff > 0
        sb(bff,j)=1./modv(blocchi(bff,j),4);
        count(blocchi(bff,j)) = count(blocchi(bff,j))+1;
        modv(blocchi(bff,j),5) = count(blocchi(bff,j));
    end
end

fover = modv(:,5)>=1; % We only store the blocks crossed by at least 1 ray.
inblocchi = modv(fover,:);
inblocchi(:,5)=find(fover);% select only the blocks with a number of rays >= 1
