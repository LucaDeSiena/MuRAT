function hitmap(Murat, x0, y0, x1, y1, Ns)
% 
% HITMAP produces horizontal and vertical cross section displaying the
%                   seismic ray count as hitmap. 3 default sections in the
%                   model are built (N-S, W-E, and horizontal). Than,
%                   starting and ending coordinates of a custom oblique
%                   cross sections can be added.
%
% Input parameters:
%    Murat:         Data structure as built by Murat
%    x0:            Longitude of the custom starting point for the oblique
%                   cross section
%    y0:            Latitude of the custom starting point for the oblique
%                   cross section
%    x1:            Longitude of the custom ending point for the oblique
%                   cross section
%    y1:            Latitude of the custom ending point for the oblique
%                   cross section
%    Ns:            Number of points along the custom oblique section
%
% Output parameters:
%    Horizontal_rays_hitmap.png:  Horizontal cross section of the total ray
%                                 count as hitmap
%    WE_rays_hitmap.png:          Vertical W-E cross section of the total
%                                 ray count as hitmap
%    SN_rays_hitmap.png:          Vertical S-N cross section of the total
%                                 ray count as hitmap
%    Section_rays_hitmap:         Oblique vertical cross section of the total
%                                 ray count as hitmap
%

% Dafault values for starting and ending points
arguments
    Murat
    x0 = NaN
    y0 = NaN
    x1 = NaN
    y1 = NaN
    Ns = 100
end

disp('Building horizontal and vertical hitmaps')

% Get the data from Murat structure
try
    rayCrossing = Murat.rays.rayCrossing';
catch
    rayCrossing = Murat.data.rayCrossing';
end
grid = Murat.input.DDcoordinates;
x = grid(:,1);
y = grid(:,2);
z = grid(:,3);
tbl = table(x,y,z,rayCrossing);

% Build default hitmaps
figure;
h = heatmap(tbl,'x','y','ColorVariable','rayCrossing','Colormap', hot);
h.NodeChildren(3).YDir='normal';
saveas(gcf,'Horizontal_rays_hitmap.png');
figure;
vx = heatmap(tbl,'x','z','ColorVariable','rayCrossing','Colormap', hot);
vx.NodeChildren(3).YDir='normal';
saveas(gcf,'WE_rays_hitmap.png');
figure;
vy = heatmap(tbl,'y','z','ColorVariable','rayCrossing','Colormap', hot);
vy.NodeChildren(3).YDir='normal';
saveas(gcf,'SN_rays_hitmap.png');

if ~isnan(x0) && ~isnan(y0) && ~isnan(x1) && ~isnan(y1)
    % Build the new oblique vertical section
    dx = x1 - x0;
    dy = y1 - y0;
    L = sqrt(dx^2 + dy^2);
    if L > 0
        s = linspace(0, L, Ns);
        zq = unique(z);
        [Xs, Zs] = meshgrid(s, zq);
        Xs = x0 + (x1 - x0) * (Xs / L);
        Ys = y0 + (y1 - y0) * (Xs - x0) / (x1 - x0);
        
        % Value interpolation and plot
        F = scatteredInterpolant(x,y,z,rayCrossing, 'natural', 'none');
        Vs = F(Xs, Ys, Zs);
        figure;
        imagesc(s, zq, Vs);
        set(gca, 'YDir', 'normal','Colormap', hot);
        xlabel('Distance along section');
        ylabel('Elevation (m)');
        colorbar;
        saveas(gcf,'Section_rays_hitmap.png');
    else
        disp("Error: section length is 0")
    end
end

