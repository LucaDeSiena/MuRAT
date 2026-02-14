function outFile = saveFigureAsImage(pathFolder, fig)
% saveFigureAsImage  Capture figure and write to PNG.
%   outFile = saveFigureAsImage(pathFolder) captures gcf and writes to
%   pathFolder.png (or pathFolder if it already has .png).
%   outFile = saveFigureAsImage(pathFolder, fig) captures specified figure
%   handle.
%
%   Returns full path of written file.

arguments
    pathFolder (1,:) char
    fig = gcf
end

% Resolve figure handle
if isempty(fig) || ~ishandle(fig)
    fig = gcf;
    fig.Color = 'white';
end

% Ensure filename ends with .png or treat as folder
[pth, name, ext] = fileparts(pathFolder);
if isempty(ext)
    % If pathFolder is an existing folder, use default name
    if isfolder(pathFolder)
        pth = pathFolder;
        name = 'figure';
    end
    ext = '.png';
end

outFile = fullfile(pth, [name, ext]);

% Create parent folder if needed
[outDir, ~, ~] = fileparts(outFile);
if ~isempty(outDir) && ~isfolder(outDir)
    mkdir(outDir);
end

% Use a fast exporter
% origRenderer = fig.Renderer;
% fig.Renderer = 'painter';        % good for raster PNGs
ax = findall(fig, 'Type', 'axes');
for k = 1:numel(ax)
    ax(k).Color = 'white';      % axes background
    ax(k).XColor = 'k';         % tick and axis line color
    ax(k).YColor = 'k';
    if isprop(ax(k), 'ZColor')
        ax(k).ZColor = 'k';
    end
    % Title and labels (if present)
    if ~isempty(ax(k).Title)
        ax(k).Title.Color = 'k';
    end
    if ~isempty(ax(k).XLabel)
        ax(k).XLabel.Color = 'k';
    end
    if ~isempty(ax(k).YLabel)
        ax(k).YLabel.Color = 'k';
    end
end

% Make all text objects black (covers text(), annotation(), etc.)
txt = findall(fig, 'Type', 'text');
set(txt, 'Color', 'k');

% Configure all legends: white box background and black text
leg = findall(fig, 'Type', 'legend');
for k = 1:numel(leg)
    % Legend text color
    leg(k).TextColor = 'k';
    % Legend background (opaque white)
    if isprop(leg(k), 'Color')
        leg(k).Color = 'white';    % older MATLAB also works
    end
    % For newer MATLAB, set the box face color to ensure opacity
    if isprop(leg(k), 'BoxFace')
        try
            leg(k).BoxFace.ColorType = 'truecoloralpha'; % ensure supports alpha
            leg(k).BoxFace.Color = [1 1 1];              % white
            leg(k).BoxFace.FaceAlpha = 1;               % fully opaque
        catch
            % ignore if property not supported
        end
    end
end

% Update all colorbars
cbs = findall(fig, 'Type', 'colorbar');
for k = 1:numel(cbs)
    cb = cbs(k);
    % Tick and label color (affects tick labels and tick marks)
    if isprop(cb, 'Color')
        cb.Color = 'k';        % black tick labels and tick marks
    else
        set(cb, 'Color', 'k');
    end

    % Label (colorbar title) – make sure it's black
    if ~isempty(cb.Label)
        cb.Label.Color = 'k';
    end

    % If legend-like box face exists (newer MATLAB), make it opaque white
    if isprop(cb, 'BoxFace')
        try
            cb.BoxFace.ColorType = 'truecoloralpha';
            cb.BoxFace.Color = [1 1 1];   % white
            cb.BoxFace.FaceAlpha = 1;     % fully opaque
        catch
            % ignore if not supported
        end
    else
        % Fallback: try setting the colorbar background
        try
            cb.Color = 'k';              % keep ticks black
            % no direct background property in older releases; use axes color
            % set parent axes background white if needed:
            ax = cb.Parent;
            if isprop(ax, 'Color')
                ax.Color = 'white';
            end
        catch
        end
    end
end


drawnow;                        % ensure content is up-to-date

% Use exportgraphics which is typically faster
% You can set resolution with 'Resolution',300 if needed.
exportgraphics(fig, outFile, BackgroundColor='white', Resolution=72);

% restore renderer
% fig.Renderer = origRenderer;

end