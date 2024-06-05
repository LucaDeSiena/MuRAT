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

locationsDeg                =   Murat.rays.locationsDeg;
uncertaintyQc               =   Murat.Qc.uncertaintyQc;
Murat.rays.locDegOriginal   =   locationsDeg;

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
Murat.Q.energyRatioBodyCoda =...
                            Murat.Q.energyRatioBodyCoda(ind_to_keep,:);
Murat.Q.energyRatioCodaNoise=...
    Murat.Q.energyRatioCodaNoise(ind_to_keep,:);
Murat.Qc.inverseQc          =   Murat.Qc.inverseQc(ind_to_keep,:);
Murat.PD.inversionMatrixPeakDelay...
                            =...
    Murat.PD.inversionMatrixPeakDelay(ind_to_keep,:);
Murat.Q.inversionMatrixQ    =   Murat.Q.inversionMatrixQ(ind_to_keep,:);
Murat.Qc.inversionMatrixQc  =   Murat.Qc.inversionMatrixQc(ind_to_keep,:);
Murat.rays.locationsDeg     =   Murat.rays.locationsDeg(ind_to_keep,:);
Murat.PD.peakDelay          =   Murat.PD.peakDelay(ind_to_keep,:);
Murat.PD.retainPeakDelay    =   Murat.PD.retainPeakDelay(ind_to_keep,:);
Murat.Q.retainQ             =   Murat.Q.retainQ(ind_to_keep,:);
Murat.Qc.retainQc           =   Murat.Qc.retainQc(ind_to_keep,:);
Murat.Qc.tCoda              =   Murat.Qc.tCoda(ind_to_keep,:);
Murat.rays.totalLengthRay   =   Murat.rays.totalLengthRay(ind_to_keep,:);
Murat.rays.travelTime       =   Murat.rays.travelTime(ind_to_keep,:);
Murat.Qc.uncertaintyQc      =   Murat.Qc.uncertaintyQc(ind_to_keep,:);
Murat.PD.variationPeakDelay =   Murat.PD.variationPeakDelay(ind_to_keep,:);

disp(['used waveforms after declustering: ',...
    num2str(length(ind_to_keep)*components)])
end