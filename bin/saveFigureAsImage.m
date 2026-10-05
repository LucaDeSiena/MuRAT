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

% Force white figure background
fig.Color = [1 1 1];

% Find all axes (including tiled, polar, etc.)
ax = findall(fig, 'Type', 'axes');
for k = 1:numel(ax)
    % Axes background white
    ax(k).Color = [1 1 1];
    % Axes lines, tick labels, and title/label text to black
    ax(k).XColor = [0 0 0];
    ax(k).YColor = [0 0 0];
    if isprop(ax(k),'ZColor')       % cartesian 3D axes
        ax(k).ZColor = [0 0 0];
    end
    % Set label and title colors (covers xlabel/ylabel/title objects)
    try
        ax(k).Title.Color = [0 0 0];
    catch
    end
    try
        ax(k).XLabel.Color = [0 0 0];
        ax(k).YLabel.Color = [0 0 0];
        if isprop(ax(k),'ZLabel')
            ax(k).ZLabel.Color = [0 0 0];
        end
    catch
    end
    % Tick label interpreter may be preserved; ensure color for text objects
    tickText = findall(ax(k), 'Type', 'text');
    for t = 1:numel(tickText)
        tickText(t).Color = [0 0 0];
    end
end

% Configure legends
lg = findall(fig, 'Type', 'legend');
for L = 1:numel(lg)
    lg(L).TextColor   = [0 0 0];   % legend text black
    lg(L).Color       = [1 1 1];   % legend background white
    lg(L).EdgeColor   = [0 0 0];   % legend border black
    % If legend contains title
    try
        lg(L).Title.Color = [0 0 0];
    catch
    end
end

drawnow; % ensure updates take effect

% Resolve figure handle
if isempty(fig) || ~ishandle(fig)
    fig = gcf;
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
origRenderer = fig.Renderer;
fig.Renderer = 'painter';        % good for raster PNGs
drawnow;                        % ensure content is up-to-date

% Use exportgraphics which is typically faster
% You can set resolution with 'Resolution',300 if needed.
exportgraphics(fig, outFile, 'BackgroundColor', 'white', 'Resolution', 72);

% restore renderer
fig.Renderer = origRenderer;

end