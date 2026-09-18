function f = myfig(figName, varargin)
% myfig  Create figure with preset defaults and optional name.
%   f = myfig()                 -> creates figure with defaults, no title.
%   f = myfig('My Title')       -> creates figure with Name = 'My Title'.
%   f = myfig('', 'Visible','on','Position',[0 0 800 600]) 
%                               -> override defaults with Name-Value pairs.
%
%   Returns the Figure object.

if nargin < 1, figName = ''; end

% default properties
defaultPos = [20 400 1200 1000];
defaults = { ...
    'NumberTitle', 'off', ...
    'Position',    defaultPos, ...
    'Visible',     'off' ...
    };

% include Name if provided
if ~isempty(figName)
    defaults = [{'Name', figName}, defaults];
end

% call figure with defaults first, then user overrides (varargin last wins)
f = figure(defaults{:}, varargin{:});

end