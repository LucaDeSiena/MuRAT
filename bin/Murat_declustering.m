function Murat = Murat_declustering(Murat,factor)

% function Murat = Murat_declustering(Murat,factor)

%PLOTS rays before and after declustering

% Input parameters:
% Murat structure variable
% factor:                   factor by which original grid is divided into

disp(['used waveforms before declustering: ', num2str(length(Murat.data.locationsDeg)*Murat.input.components)])

% plot events before declustering to see what is there
figure
hold on
for n = 1:length(Murat.data.locationsDeg)
    plot([Murat.data.locationsDeg(n,2),Murat.data.locationsDeg(n,5)],...
        [Murat.data.locationsDeg(n,1),Murat.data.locationsDeg(n,4)],'-k')
end
title('comparison between original and declustered data set')
ylim([Murat.input.origin(1),Murat.input.end(1)])
xlim([Murat.input.origin(2),Murat.input.end(2)])
xlabel('longitude [°]')
ylabel('latitude [°]')
grid on

% create smaller grid to cluster events per grid
lat_step = length(Murat.input.y)*factor;
lon_step = length(Murat.input.x)*factor;
z_step   = length(Murat.input.z)*factor;
lat_grid = linspace(Murat.input.origin(2),Murat.input.end(2),lat_step);
lon_grid = linspace(Murat.input.origin(1),Murat.input.end(1),lon_step);
z_grid   = -linspace(Murat.input.origin(3),Murat.input.end(3),z_step);

% decimate evestaz to one entry per event - station pair
evestaz = unique(Murat.data.locationsDeg,'rows','stable');

%% loop around new grid
new_evestaz = [];

for i=1:lon_step-1
    
    for j = 1:lat_step-1
        
        for k = 1:z_step-1
            % find all events within one grid cell
            find_evs = find(evestaz(:,1)>lon_grid(i) & ...
                evestaz(:,1)<lon_grid(i+1) & ...
                evestaz(:,2)>lat_grid(j) & ...
                evestaz(:,2)<lat_grid(j+1) & ...
                evestaz(:,3)>z_grid(k) & ...
                evestaz(:,3)<z_grid(k+1));
            
            % check out each grid cell and only keep events with highest RZZ
            if ~isempty(find_evs)
                % get RZZ & events/stations & indice
                events = evestaz(find_evs,:);
                events(:,7) = Murat.data.uncertaintyQc(find_evs);
                events(:,8) = find_evs;
                % check if stations are double, if so, only keep event/station
                % pair with highest RZZ
                d = sortrows(events, [4 7]);
                [~, ia, ~] = unique(d(:,4),'rows','last');
                to_keep = d(ia,:);
                new_evestaz = [new_evestaz;to_keep];
            end
        end
    end
    
end

new_evestaz = sortrows(new_evestaz,8);

%% remove deleted data from Murat structure variable
ind_to_keep                     =   new_evestaz(:,8);
Murat.data.energyRatioBodyCoda  =...
    Murat.data.energyRatioBodyCoda(ind_to_keep,:);
Murat.data.energyRatioCodaNoise =...
    Murat.data.energyRatioCodaNoise(ind_to_keep,:);
Murat.data.inverseQc            =   Murat.data.inverseQc(ind_to_keep,:);
Murat.data.inversionMatrixPeakDelay      =...
    Murat.data.inversionMatrixPeakDelay(ind_to_keep,:);
Murat.data.inversionMatrixQ     =...
    Murat.data.inversionMatrixQ(ind_to_keep,:);
Murat.data.inversionMatrixQc    =...
    Murat.data.inversionMatrixQc(ind_to_keep,:);
Murat.data.locationsDeg         =   Murat.data.locationsDeg(ind_to_keep,:);
Murat.data.peakDelay            =   Murat.data.peakDelay(ind_to_keep,:);
Murat.data.retainPeakDelay      =   Murat.data.retainPeakDelay(ind_to_keep,:);
Murat.data.retainQ              =   Murat.data.retainQ(ind_to_keep,:);
Murat.data.retainQc             =   Murat.data.retainQc(ind_to_keep,:);
Murat.data.tCoda                =   Murat.data.tCoda(ind_to_keep,:);
Murat.data.totalLengthRay       =   Murat.data.totalLengthRay(ind_to_keep,:);
Murat.data.travelTime           =   Murat.data.travelTime(ind_to_keep,:);
Murat.data.uncertaintyQc        =   Murat.data.uncertaintyQc(ind_to_keep,:);
Murat.data.variationPeakDelay   =   Murat.data.variationPeakDelay(ind_to_keep,:);
% triple entries for evestaz to match original list
% k = 1;
% Murat.data.locationsDeg = [];
% for n = 1:length(new_evestaz)
%     Murat.data.locationsDeg(k,:)=   new_evestaz(n,1:6);
%     k = k+1;
%     Murat.data.locationsDeg(k,:)=   new_evestaz(n,1:6);
%     k = k+1;
%     Murat.data.locationsDeg(k,:)=   new_evestaz(n,1:6);
%     k = k+1;
% end

% check which events were kept before and compare with new list
% for i=1:3
%     find_common                     =   ismember(retainQm(:,i),ind_to_keep);
%     Murat.data.retainQc(:,i)        =   Murat.data.retainQc(find_common,i);
% end
% % get new indice
% retainQM = Murat.data.retainQc; Murat.data.retainQm = [];
% for n = 1:length(retainQM)
%     new_index(n) = find(ind_to_keep == retainQM(n));
% end
% new_index = new_index';
% Murat.data.retainQm = new_index;
% % check which events were discarded before and compare with new list
% find_common = ismember(Murat.data.discardQm,ind_to_keep);
% Murat.data.discardQm            =   Murat.data.discardQm(find_common);
% % get new indice
% discardQM = Murat.data.discardQm; Murat.data.discardQm = [];
% for n = 1:length(discardQM)
%     new_index2(n) = find(ind_to_keep == discardQM(n));
% end
% new_index2 = new_index2';
% Murat.data.discardQm = new_index2;

% Murat.data.fitrobust          =   % keep or update?
% Murat.data.averageQc          =    % keep or update?

%% plot again to check
for n = 1:length(Murat.data.locationsDeg)
    plot([Murat.data.locationsDeg(n,2),Murat.data.locationsDeg(n,5)],...
        [Murat.data.locationsDeg(n,1),Murat.data.locationsDeg(n,4)],'-r')
end
hold off

%%
disp(['used waveforms after declustering: ', num2str(length(ind_to_keep)*Murat.input.components)])
end