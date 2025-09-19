function plot_illumination(T, FF, step_deg, parallel, sunFlag)
% plot_illumination plots terminator and eclipse shadows using unified dataTable approach
%
%   Inputs:
%       T (double): time in mjd2000 days.
%       FF (double): Parameter for eclipse calculation
%       step_deg (double, optional): Step size for plotting. Default is 2.
%       parallel (logical, optional): Flag for parallel processing. Default is false.
%       sunFlag (logical, optional): Flag to plot the subsolar point. Default is false.

% Set default values
if nargin < 5, sunFlag = false; end
if nargin < 4, parallel = false; end
if nargin < 3, step_deg = 2; end
if nargin < 2, FF = 0; end

% Initialize arrays to hold graphics handles and their labels for the legend
legendHandles = [];
legendLabels = {};

% Convert MJD2000 to a normal date string
dateVec = mjd20002date(T);
normalDateStr = sprintf('%04d-%02d-%02d %02d:%02d:%02.0f', ...
    dateVec(1), dateVec(2), dateVec(3), ...
    dateVec(4), dateVec(5), dateVec(6));

% Draw the base Mercator texture
plotTextureLatLong(1, 6);
hold on;

%% Create the grid and dataTable here
step = deg2rad(step_deg);
lat_vector = (-pi/2:step:pi/2)';
lon_vector = (0:step:2*pi)';
[lat_grid, lon_grid] = meshgrid(lat_vector, lon_vector);
T_array = T * ones(numel(lat_grid), 1);

% Create dataTable with the grid
dataTable = table(lat_grid(:), lon_grid(:), T_array, ...
    'VariableNames', {'lat', 'lon', 'T'});

%% Call what_illumination with the dataTable
dataTable = what_illumination(dataTable, FF, parallel);

%% Plot based on dataTable
if FF < 10000
    % For FF < 10000, plot both terminator and eclipse
    
    % Plot nightside (terminator)
    night_indices = find(strcmp(dataTable.terminator, 'night'));
    if ~isempty(night_indices)
        hNightside = plotNightsideFromData(dataTable, night_indices);
        legendHandles(end+1) = hNightside;
        legendLabels{end+1} = 'Nightside';
    end
    
    % Plot penumbra
    penumbra_indices = find(strcmp(dataTable.eclipse, 'Penumbra'));
    if ~isempty(penumbra_indices)
        hPenumbra = plotEclipseRegionFromData(dataTable, penumbra_indices, 'Penumbra');
        if ~isempty(hPenumbra)
            legendHandles(end+1) = hPenumbra;
            legendLabels{end+1} = 'Penumbra';
        end
    end
    
    % Plot umbra
    umbra_indices = find(strcmp(dataTable.eclipse, 'Umbra'));
    if ~isempty(umbra_indices)
        hUmbra = plotEclipseRegionFromData(dataTable, umbra_indices, 'Umbra');
        if ~isempty(hUmbra)
            legendHandles(end+1) = hUmbra;
            legendLabels{end+1} = 'Umbra';
        end
    end
    
    title(sprintf('Terminator and Eclipse Shadows on %s', normalDateStr));
else
    % For FF >= 10000, plot only the terminator
    
    % Plot nightside (terminator)
    night_indices = find(strcmp(dataTable.terminator, 'night'));
    if ~isempty(night_indices)
        hNightside = plotNightsideFromData(dataTable, night_indices);
        legendHandles(end+1) = hNightside;
        legendLabels{end+1} = 'Nightside';
    end
    
    title(sprintf('Terminator Shadow on %s', normalDateStr));
end

% Plot sub-solar point if requested
if sunFlag
    [lat_sun, lon_sun] = sunSubPointOnEnceladus(T);
    lon_sun = mod(lon_sun + 360, 360);
    hSun = plot(lon_sun, lat_sun, 'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
    legendHandles(end+1) = hSun;
    legendLabels{end+1} = 'Sub-solar point';
end

% Create the final, consolidated legend
if ~isempty(legendHandles)
    legend(legendHandles, legendLabels, 'Location', 'best');
end

xlim([0 360]);
ylim([-90 90]);
hold off;

end

%% Helper function to plot nightside from dataTable
function hPatch = plotNightsideFromData(dataTable, night_indices)
    % Get the terminator boundary
    lat_gt = dataTable.lat(night_indices);
    lon_gt = dataTable.lon(night_indices);
    
    % Convert to degrees
    lat_deg = rad2deg(lat_gt);
    lon_deg = rad2deg(lon_gt);
    
    % Create nightside patch
    % For simplicity, we'll create a patch that covers the night region
    % This is a simplified approach - you may need to refine based on your exact needs
    [lat_sun, ~] = sunSubPointOnEnceladus(dataTable.T(1));
    
    if lat_sun >= 0
        edgeLat = -90;
    else
        edgeLat = +90;
    end
    
    % Create terminator curve (simplified approach)
    npts = 361;
    lon_plot = linspace(0, 360, npts);
    [~, lon_sun] = sunSubPointOnEnceladus(dataTable.T(1));
    lon_sun = mod(lon_sun + 360, 360);
    lha = deg2rad(lon_plot - lon_sun);
    lat_formula = atan(-cos(lha) ./ tan(deg2rad(lat_sun)));
    lat = rad2deg(lat_formula);
    
    px = [lon_plot, fliplr(lon_plot)];
    py = [lat, edgeLat*ones(1, npts)];
    hPatch = patch(px, py, 'black', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
end

%% Helper function to plot eclipse regions from dataTable
function hPatch = plotEclipseRegionFromData(dataTable, indices, type)
    hPatch = [];
    if isempty(indices)
        return;
    end
    
    % Set colors based on type
    if strcmp(type, 'Penumbra')
        color = [0.6 0.4 0.2];
        alpha = 0.4;
    else % Umbra
        color = [0.4 0.2 0];
        alpha = 0.6;
    end
    
    % Process first longitude part
    indices_part1 = indices(dataTable.lon(indices) <= pi/2);
    if ~isempty(indices_part1)
        lat_deg = rad2deg(dataTable.lat(indices_part1));
        lon_deg = rad2deg(dataTable.lon(indices_part1));
        if length(unique(lon_deg)) > 1 && length(unique(lat_deg)) > 1
            k = boundary(lon_deg, lat_deg, 1);
            p = patch(lon_deg(k), lat_deg(k), color, 'FaceAlpha', alpha, 'EdgeColor', 'none');
            if isempty(hPatch), hPatch = p; end
        end
    end
    
    % Process second longitude part
    indices_part2 = indices(dataTable.lon(indices) >= 3/2*pi);
    if ~isempty(indices_part2)
        lat_deg = rad2deg(dataTable.lat(indices_part2));
        lon_deg = rad2deg(dataTable.lon(indices_part2));
        if length(unique(lon_deg)) > 1 && length(unique(lat_deg)) > 1
            k = boundary(lon_deg, lat_deg, 1);
            p = patch(lon_deg(k), lat_deg(k), color, 'FaceAlpha', alpha, 'EdgeColor', 'none');
            if isempty(hPatch), hPatch = p; end
        end
    end
end