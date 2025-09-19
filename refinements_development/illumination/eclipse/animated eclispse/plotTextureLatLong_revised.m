function [fig] = plotTextureLatLong_revised(bodyName, ax)
% This function plots a body's texture map onto a specified axes.
% It is revised for use in animations to prevent creating new figures.

% Author: Andrea Bellome, Revised for Animation
% Last revision: 28/08/2025

%% INPUTS %%
% - bodyName : The name of the body as a string (e.g., 'Enceladus').
%              This name should match the texture file (e.g., 'Enceladus.jpg').
% - ax       : (Optional) Handle to the axes on which to plot. If not
%              provided, a new figure and axes will be created.

%% OUTPUTS %%
% - fig : The figure handle.

%% FUNCTION %%

% If no axes handle is provided, create a new figure and axes.
if nargin < 2 || isempty(ax)
    fig = figure('Color', [1 1 1]);
    ax = axes('Parent', fig);
else
    fig = get(ax, 'Parent'); % Get the figure handle from the axes
end

% Set the provided axes as the current axes for plotting
axes(ax);
hold on; 

% Set up the axes properties
xlabel(ax, 'Longitude - deg');
ylabel(ax, 'Latitude - deg');
axis(ax, 'equal');
xlim(ax, [0 360]); 
ylim(ax, [-90 90]);
xticks(ax, 0:60:360);
xticklabels(ax, {'0°','60°','120°','180°','240°','300°','360°'});
grid(ax, 'on');

try
    textureFile = [bodyName '.jpg'];
    I = imread(textureFile);
catch ME
    error('Could not read texture file "%s". Make sure it is in the MATLAB path. Original error: %s', textureFile, ME.message);
end

% Use the image function on the specified axes
h = image(ax, xlim(ax), -ylim(ax), I);

% Push the texture map to the background
uistack(h, 'bottom');

hold(ax, 'off');

end