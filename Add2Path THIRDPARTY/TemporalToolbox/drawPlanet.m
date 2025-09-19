function [hfig] = drawPlanet( namePlanet, position, handle, scale)
% drawPlanet.m - draws a planet surface (approx) in a given position.
%
% PROTOTYPE: 
%           [HFIG] = drawPlanet( namePlanet, position, varargin)
%
% DESCRIPTION: 
%  The function draws a planet in a given position and also a
%  given axis.
%
% INPUT:
%  namePlanet     -    string with the name of the planet (not case
%                           sensitive)
%  position       -    centre of the planet to draw
%  handle         -    Figurr or axis handle

%
% OUTPUT:
%  HFIG           -    Figure object of the planet draw.
%
% EXAMPLE: HFIG=drawPlanet('Sun',[0 0 0]);
%
% CALLED FUNCTIONS: getAstroConstants.m
%
% AUTHOR:
%   Joan Pau Sanchez, 14/10/2009, MATLAB, drawPlanet.m
%   
% CHANGELOG:
% REVISION:
%--------------------------------------------------------------------------

% 1. Searching Texture Directory
persistent pathTextures
if isempty(pathTextures)
    p = path;
    tokens = strsplit(p, pathsep);    % pathsep = ':' on macOS, ';' on Windows
    pathTextures = '';
    for k = 1:numel(tokens)
        if contains(tokens{k}, 'textures')
            pathTextures = tokens{k};
            break
        end
    end
    if isempty(pathTextures)
        error('TOOLBOX:%s:propertyError','textures folder not on path');
    end
end
%--------------------------------------------------------------------------
% 2. Loading Texture Map 
     files = dir(pathTextures);
     found = false;
for i=1:numel(files)
    if contains(files(i).name, namePlanet)
        [texturemap,map] = imread( fullfile(pathTextures, files(i).name) );
        found = true;
        break
    end
end
if ~found
    error('TOOLBOX:%s:propertyError','couldn''t find texture for %s',namePlanet);
end
%--------------------------------------------------------------------------
% 3. Preraring Figure Object
if nargin<3
    HAXIS = gca;
elseif ishandle(handle)==0
        msg = ['The figure handle is not valid'];
        eid = sprintf('TOOLBOX:%s:propertyError', mfilename);
        error(eid,'%s',msg)
else
    try
        HAXIS=gca(handle);
    catch
        HAXIS=handle;  
    end
    hold on
end
%--------------------------------------------------------------------------
if nargin<4
    scale=1;   
end
% 4. Ploting Planet
[radPlanet] = getAstroConstants(namePlanet, 'Radius');
     nfaces = 50;
    [X,Y,Z] = sphere(nfaces);
% 4.1. scaling the sphere and locating the planet
X = -radPlanet*X*scale + position(1);
Y = -radPlanet*Y*scale + position(2);
Z = -radPlanet*Z*scale + position(3);
% 4.2. ploting
surf(HAXIS, X, Y, Z,texturemap, 'LineStyle', 'none'...
                                       , 'FaceColor', 'texturemap');
hfig = gcf; 
return