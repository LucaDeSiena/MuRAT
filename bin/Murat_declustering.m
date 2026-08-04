function indKeep =  Murat_declustering(origin,ending,x,y,z,QcM,uncertainty,locationsDeg,factor)
% indKeep =  Murat_declustering(origin,ending,x,y,z,QcM,uncertainty,locationsDeg,factor)
%
% SAVES declustered variables efore and after declustering.
% The declustering is based on dividing the original grid into finer node
% spacing and selecting those with the least coda uncertainty. 
%
% Input parameters:
%    Murat:                Murat structure variable
%    factor:               input factor used to divide the original grid

% create smaller grid to cluster events
lat_step    =	length(y)*factor;
lon_step    =   length(x)*factor;
z_step      =   length(z)*factor;

lat_grid    =   linspace(origin(2),ending(2),lat_step);
lon_grid    =   linspace(origin(1),ending(1),lon_step);
z_grid      =   -linspace(origin(3),ending(3),z_step);

% decimate locationsDeg to one entry per event - station pair
locationsDeg    =   unique(locationsDeg,'rows','stable');

indKeep         =   false(length(uncertainty),1);
new_locationsDeg=   [];

for i = 1:lon_step-1
    
    for j = 1:lat_step-1
        
        for k = 1:z_step-1
            % find all events within one grid cell
            find_evs    =   find(locationsDeg(:,1)>lon_grid(i) & ...
                locationsDeg(:,1)<lon_grid(i+1) & ...
                locationsDeg(:,2)>lat_grid(j) &...
                locationsDeg(:,2)<lat_grid(j+1) & ...
                locationsDeg(:,3)>z_grid(k) &...
                locationsDeg(:,3)<z_grid(k+1));
            
            % check out each grid cell, only keep events with highest RZZ
            if ~isempty(find_evs)
                % get RZZ & events/stations & indice
                events      =   locationsDeg(find_evs,:);

                if isequal(QcM,'Linear')
                    events(:,7) =   uncertainty(find_evs);
                else
                    events(:,7) =   1./uncertainty(find_evs);
                end
                events(:,8)     =	find_evs;
                % check if stations are double, if so, only keep
                % event/station pair with lowest uncertainty on Qc
                d               =   sortrows(events, [4 7]);
                [~, ia, ~]      =   unique(d(:,4),'rows','last');
                to_keep         =   d(ia,:);
                new_locationsDeg=   [new_locationsDeg;to_keep]; %#ok<AGROW>
            end
        end
    end
    
end

new_locationsDeg                =   sortrows(new_locationsDeg,8);

indKeep(new_locationsDeg(:,8))  =   true;

end