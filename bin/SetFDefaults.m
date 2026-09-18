function SetFDefaults()
% DEFAULT settings for MuRAT figures
h = gca;
if isa(h,'matlab.graphics.axis.Axes')
    safeSet(h, 'GridColor', [0 0 0]);
    safeSet(h, 'GridLineStyle', '--');
    safeSet(h, 'GridAlpha', 0.3);
    safeSet(h, 'LineWidth', 1.5);
    safeSet(h, 'FontSize', 12);
    safeSet(h, 'XGrid', 'on');
    safeSet(h, 'YGrid', 'on');
end
grid on

% Helper: set property only when supported and catch set failures
function safeSet(h, prop, val)
if any(strcmp(prop, properties(h)))
    try
        set(h, prop, val);
    catch
        % ignore failures silently
    end
end