function Murat = Murat_declustering(Murat,factor)
% function Murat = Murat_declustering(Murat,factor)
%
% PLOTS rays before and after declustering
%
% Input parameters:
%    Murat:                Murat structure variable
%    factor:               factor by which original grid is divided into
components                  =   Murat.input.components;
origin                      =   Murat.input.origin;
ending                      =   Murat.input.end;
x                           =   Murat.input.x;
y                           =   Murat.input.y;
z                           =   Murat.input.z;

locationsDeg                =   Murat.data.locationsDeg;
uncertaintyQc               =   Murat.data.uncertaintyQc;
Murat.data.locDegOriginal   =   locationsDeg;

disp(['used waveforms before declustering: ',...
    num2str(length(locationsDeg)*components)])

% create smaller grid to cluster events
lat_step                    =	length(y)*factor;
lon_step                    =   length(x)*factor;
z_step                      =   length(z)*factor;
lat_grid                    =   linspace(origin(2),ending(2),lat_step);
lon_grid                    =   linspace(origin(1),ending(1),lon_step);
z_grid                      =   -linspace(origin(3),ending(3),z_step);

% decimate evestaz to one entry per event - station pair
evestaz                     =   unique(locationsDeg,'rows','stable');

%% loop around new grid
new_evestaz                 =   [];

for i = 1:lon_step-1
    
    for j = 1:lat_step-1
        
        for k = 1:z_step-1
            % find all events within one grid cell
            find_evs        =   find(evestaz(:,1)>lon_grid(i) & ...
                evestaz(:,1)<lon_grid(i+1) & ...
                evestaz(:,2)>lat_grid(j) & evestaz(:,2)<lat_grid(j+1) & ...
                evestaz(:,3)>z_grid(k) & evestaz(:,3)<z_grid(k+1));
            
            % check out each grid cell, only keep events with highest RZZ
            if ~isempty(find_evs)
                % get RZZ & events/stations & indice
                events      =   evestaz(find_evs,:);
                events(:,7) =   uncertaintyQc(find_evs);
                events(:,8) =	find_evs;
                % check if stations are double, if so, only keep
                % event/station pair with highest RZZ
                d           =   sortrows(events, [4 7]);
                [~, ia, ~]  =   unique(d(:,4),'rows','last');
                to_keep     =   d(ia,:);
                new_evestaz =   [new_evestaz;to_keep]; %#ok<AGROW>
            end
        end
    end
    
end

new_evestaz                 =   sortrows(new_evestaz,8);

%% remove deleted data from Murat structure variable
ind_to_keep                 =   new_evestaz(:,8);
Murat.data.energyRatioBodyCoda  =...
    Murat.data.energyRatioBodyCoda(ind_to_keep,:);
Murat.data.energyRatioCodaNoise =...
    Murat.data.energyRatioCodaNoise(ind_to_keep,:);
Murat.data.inverseQc        =   Murat.data.inverseQc(ind_to_keep,:);
Murat.data.inversionMatrixPeakDelay	=...
    Murat.data.inversionMatrixPeakDelay(ind_to_keep,:);
Murat.data.inversionMatrixQ =...
    Murat.data.inversionMatrixQ(ind_to_keep,:);
Murat.data.inversionMatrixQc=...
    Murat.data.inversionMatrixQc(ind_to_keep,:);
Murat.data.locationsDeg     =   Murat.data.locationsDeg(ind_to_keep,:);
Murat.data.peakDelay        =   Murat.data.peakDelay(ind_to_keep,:);
Murat.data.retainPeakDelay  =   Murat.data.retainPeakDelay(ind_to_keep,:);
Murat.data.retainQ          =   Murat.data.retainQ(ind_to_keep,:);
Murat.data.retainQc         =   Murat.data.retainQc(ind_to_keep,:);
Murat.data.tCoda            =   Murat.data.tCoda(ind_to_keep,:);
Murat.data.totalLengthRay   =   Murat.data.totalLengthRay(ind_to_keep,:);
Murat.data.travelTime       =   Murat.data.travelTime(ind_to_keep,:);
Murat.data.uncertaintyQc    =   Murat.data.uncertaintyQc(ind_to_keep,:);
Murat.data.variationPeakDelay =...
    Murat.data.variationPeakDelay(ind_to_keep,:);

disp(['used waveforms after declustering: ',...
    num2str(length(ind_to_keep)*components)])
end