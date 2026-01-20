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