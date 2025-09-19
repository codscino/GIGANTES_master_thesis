function [] = Plot_Flyby_GT_illuminated(Flyby_Data, color, T0, pars, flyby_duration_minutes, ...
    use_parallel, coloured_illumination, background_illumination, step_deg, incidence_angle_ranges)
% Plots the Ground-track of a Flyby with optional illumination-based coloring

%% Set defaults
if nargin < 6, use_parallel = false; end
if nargin < 7, coloured_illumination = false; end
if nargin < 8, background_illumination = false; end
if nargin < 9, step_deg = 1; end
if nargin < 10, incidence_angle_ranges = []; end

%% Extract Flyby Data
lats = Flyby_Data.lats;
longs = Flyby_Data.longs;
rp_lat = Flyby_Data.rp_lat;
rp_long = Flyby_Data.rp_long;
lats_deg = rad2deg(lats);
longs_deg = rad2deg(longs);
n_points = length(lats);

%% Calculate time array for trajectory
T_array = calculateTrajectoryTimes(T0, n_points, flyby_duration_minutes);

%% Force colored illumination if background requested
if background_illumination
    coloured_illumination = true;
end

%% Calculate illumination conditions if needed
if coloured_illumination || background_illumination
    % Create dataTable for trajectory points
    dataTable = table(lats(:), longs(:), T_array(:), 'VariableNames', {'lat', 'lon', 'T'});
    
    % Check if eclipse is possible
    is_eclipse_possible = checkEclipsePossibility(T0);
    
    % Calculate illumination
    dataTable_illum = what_illumination(dataTable, is_eclipse_possible, use_parallel);
    
    % Determine colors
    colors_array = determineColors(dataTable_illum, color, n_points, ...
        incidence_angle_ranges, lats, longs, T0);
else
    colors_array = repmat(color, n_points, 1);
end

%% Create plots
if background_illumination
    % First subplot: Colored groundtrack
    subplot(2, 1, 1);
    setupSubplot(T0, 'Flyby Groundtrack with Illumination-Based Coloring');
    plotGroundtrack(longs_deg, lats_deg, longs, colors_array, color, ...
        rp_long, rp_lat, coloured_illumination);
    xlim([0 360]);  % Ensure proper limits for this subplot
    ylim([-90 90]);
    
    % Second subplot: Illumination overlay
    subplot(2, 1, 2);
    setupSubplot(T0, 'Groundtrack with Illumination Context Overlay');
    plotIlluminationOverlay(T0, step_deg, is_eclipse_possible);
    
    % Add uniform-colored groundtrack
    colors_uniform = repmat(color, n_points, 1);
    plotGroundtrack(longs_deg, lats_deg, longs, colors_uniform, color, ...
        rp_long, rp_lat, false);
    xlim([0 360]);  % Ensure proper limits for this subplot
    ylim([-90 90]);
else
    % Single plot (no subplots)
    setupSubplot(T0, 'Flyby Groundtrack');
    plotGroundtrack(longs_deg, lats_deg, longs, colors_array, color, ...
        rp_long, rp_lat, coloured_illumination);
    xlim([0 360]);
    ylim([-90 90]);
end

end

%% === HELPER FUNCTIONS === %%

function T_array = calculateTrajectoryTimes(T0, n_points, flyby_duration_minutes)
    % SIMPLIFIED: Assume constant velocity flyby (nearly straight line trajectory)
    
    % Convert duration to days
    half_duration_days = (flyby_duration_minutes / 2) / (24 * 60);
    
    % Simply distribute time points uniformly
    % T0 is at pericentre (middle of the flyby)
    T_array = linspace(T0 - half_duration_days, T0 + half_duration_days, n_points)';
end

function is_eclipse_possible = checkEclipsePossibility(T0)
    % Check if eclipses are possible based on:
    % 1. Proximity to Saturn's equinox (within ±3 years)
    % 2. Sun-Saturn-Enceladus phase angle (165° to 195°)
    
    %% First check: Proximity to equinox
    % Saturn equinox data (in mjd2000 days)
    spring_equinox_ref = 3510;  % Reference spring equinox
    autumn_equinox_ref = 9257;  % Reference autumn equinox
    equinox_period = 5747;      % Days between equinoxes (~15.7 years)
    
    % Window around equinox where eclipses can occur
    equinox_window = 3 * 365.25;  % ±3 years in days
    
    % Find the nearest equinox to T0
    is_near_equinox = checkNearEquinox(T0, spring_equinox_ref, autumn_equinox_ref, ...
                                       equinox_period, equinox_window);
    
    if ~is_near_equinox
        fprintf('Eclipse not possible: T0 is not within ±3 years of an equinox.\n');
        fprintf('  T0 = %.1f days (mjd2000)\n', T0);
        
        % Find and report the nearest equinox for user information
        [nearest_equinox, equinox_type, days_from_equinox] = findNearestEquinox(T0, ...
            spring_equinox_ref, autumn_equinox_ref, equinox_period);
        
        fprintf('  Nearest %s equinox: %.1f days (mjd2000)\n', equinox_type, nearest_equinox);
        fprintf('  Days from equinox: %.1f (%.2f years)\n', ...
                abs(days_from_equinox), abs(days_from_equinox)/365.25);
        fprintf('  Required: within %.1f days (%.1f years)\n', ...
                equinox_window, equinox_window/365.25);
        
        is_eclipse_possible = false;
        return;
    end
    
    %% Second check: Phase angle
    spiceParam.frame = 'J2000';
    spiceParam.abcorr = 'NONE';
    spiceParam.observer = '0';  % Solar System Barycenter
    
    % Get positions
    r_saturn = EphSS_car_spice2(699, T0, true, spiceParam);
    r_enceladus = EphSS_car_spice2(602, T0, true, spiceParam);
    r_sun = EphSS_car_spice2(10, T0, true, spiceParam);
    
    % Calculate phase angle
    vec_to_enc = r_enceladus - r_saturn;
    vec_to_sun = r_sun - r_saturn;
    
    phase_angle_deg = rad2deg(acos(dot(vec_to_enc, vec_to_sun) / ...
        (norm(vec_to_enc) * norm(vec_to_sun))));
    
    % Check if phase angle allows eclipse
    is_phase_ok = (phase_angle_deg >= 165) && (phase_angle_deg <= 195);
    
    %% Combined result
    is_eclipse_possible = is_near_equinox && is_phase_ok;
    
    % Report results
    if is_eclipse_possible
        fprintf('Eclipse IS possible:\n');
        fprintf('  Near equinox: YES (within ±3 years)\n');
        fprintf('  Phase angle: %.2f° (within 165°-195° range)\n', phase_angle_deg);
    else
        fprintf('Eclipse NOT possible:\n');
        fprintf('  Near equinox: %s\n', iff(is_near_equinox, 'YES', 'NO'));
        fprintf('  Phase angle: %.2f° (%s 165°-195° range)\n', ...
                phase_angle_deg, iff(is_phase_ok, 'within', 'outside'));
    end
end

function is_near = checkNearEquinox(T0, spring_ref, autumn_ref, period, window)
    % Check if T0 is within 'window' days of any equinox
    
    % Generate list of all equinoxes within reasonable range of T0
    % We'll check ±2 periods around T0 to be safe
    n_periods_check = ceil(T0 / period) + 2;
    
    is_near = false;
    
    % Check spring equinoxes
    for n = -2:n_periods_check
        spring_equinox = spring_ref + n * period * 2;  % *2 because spring recurs every 2 periods
        if abs(T0 - spring_equinox) <= window
            is_near = true;
            return;
        end
    end
    
    % Check autumn equinoxes
    for n = -2:n_periods_check
        autumn_equinox = autumn_ref + n * period * 2;  % *2 because autumn recurs every 2 periods
        if abs(T0 - autumn_equinox) <= window
            is_near = true;
            return;
        end
    end
end

function [nearest_eq, eq_type, days_diff] = findNearestEquinox(T0, spring_ref, autumn_ref, period)
    % Find the nearest equinox to T0 and return its details
    
    % Calculate the cycle number we're in
    saturn_year = period * 2;  % Full Saturn year
    
    % Find nearest spring equinox
    n_spring = round((T0 - spring_ref) / saturn_year);
    nearest_spring = spring_ref + n_spring * saturn_year;
    dist_spring = abs(T0 - nearest_spring);
    
    % Find nearest autumn equinox  
    n_autumn = round((T0 - autumn_ref) / saturn_year);
    nearest_autumn = autumn_ref + n_autumn * saturn_year;
    dist_autumn = abs(T0 - nearest_autumn);
    
    % Return the nearest one
    if dist_spring < dist_autumn
        nearest_eq = nearest_spring;
        eq_type = 'spring';
        days_diff = T0 - nearest_spring;
    else
        nearest_eq = nearest_autumn;
        eq_type = 'autumn';
        days_diff = T0 - nearest_autumn;
    end
end

function colors_array = determineColors(dataTable_illum, base_color, n_points, ...
    incidence_ranges, lats, longs, T0)
    % Determine colors based on illumination and optional incidence angles
    
    colors_array = zeros(n_points, 3);
    
    % Apply illumination-based coloring
    for i = 1:n_points
        terminator = dataTable_illum.terminator{i};
        eclipse = dataTable_illum.eclipse{i};
        
        if strcmp(terminator, 'night')
            colors_array(i, :) = [0 0 0];  % Black
        elseif strcmp(eclipse, 'Umbra')
            colors_array(i, :) = [0.4 0.2 0];  % Dark brown
        elseif strcmp(eclipse, 'Penumbra')
            colors_array(i, :) = [0.6 0.4 0.2];  % Light brown
        else
            colors_array(i, :) = base_color;  % Original color
        end
    end
    
    % Override with incidence angle coloring if specified
    if ~isempty(incidence_ranges)
        incidence_angles = calculateIncidenceAngles(lats, longs, T0);
        outside_ranges = checkIncidenceRanges(incidence_angles, incidence_ranges);
        
        for i = 1:n_points
            if strcmp(dataTable_illum.terminator{i}, 'day') && ...
               strcmp(dataTable_illum.eclipse{i}, 'no eclipse') && ...
               outside_ranges(i) && incidence_angles(i) < 90
                colors_array(i, :) = [1 0 0];  % Red
            end
        end
    end
end

function outside = checkIncidenceRanges(angles, ranges)
    % Check which angles are outside specified ranges
    n = length(angles);
    outside = true(n, 1);
    
    for j = 1:size(ranges, 1)
        inside = (angles >= ranges(j, 1)) & (angles <= ranges(j, 2));
        outside = outside & ~inside;
    end
end

function setupSubplot(T0, titleText)
    % Setup subplot with texture and title
    plotTextureLatLong(1, 6, 1);
    
    axis equal;  % Equal scaling
    axis([0 360 -90 90]);  % Set explicit limits
    
    grid on;
    hold on;
    
    dateVec = mjd20002date(T0);
    dateStr = sprintf('%04d-%02d-%02d %02d:%02d:%02.0f', dateVec(1:6));
    title(sprintf('%s\nPericentre: %s', titleText, dateStr));
    xlabel('Longitude (deg)');
    ylabel('Latitude (deg)');
end

function plotIlluminationOverlay(T0, step_deg, is_eclipse_possible)
    % Plot nightside and eclipse regions
    
    % Plot nightside
    [lat_sun, lon_sun] = sunSubPointOnEnceladus(T0);
    lon_sun = mod(lon_sun + 360, 360);
    plotNightside(lat_sun, lon_sun);
    
    % Plot eclipse regions if possible
    if is_eclipse_possible
        plotEclipseGrid(T0, step_deg);
    end
    
    % Plot sub-solar point
    plot(lon_sun, lat_sun, 'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
end

function plotNightside(lat_sun, lon_sun)
    % Plot nightside region
    npts = 361;
    lon_plot = linspace(0, 360, npts);
    lha = deg2rad(lon_plot - lon_sun);
    lat_term = rad2deg(atan(-cos(lha) ./ tan(deg2rad(lat_sun))));
    
    edgeLat = iff(lat_sun >= 0, -90, 90);
    
    px = [lon_plot, fliplr(lon_plot)];
    py = [lat_term, edgeLat*ones(1, npts)];
    patch(px, py, 'black', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
end

function plotEclipseGrid(T0, step_deg)
    % Plot eclipse regions on a grid
    lon_range = [0:step_deg:90, 270:step_deg:360];
    lat_range = -90:step_deg:90;
    
    for lon = lon_range
        for lat = lat_range
            dataTable = table(deg2rad(lat), deg2rad(lon), T0, ...
                'VariableNames', {'lat', 'lon', 'T'});
            dataTable_illum = what_illumination(dataTable, true, false);
            
            eclipse_type = dataTable_illum.eclipse{1};
            if strcmp(eclipse_type, 'Penumbra')
                plotGridCell(lon, lat, step_deg, [0.6 0.4 0.2], 0.3);
            elseif strcmp(eclipse_type, 'Umbra')
                plotGridCell(lon, lat, step_deg, [0.4 0.2 0], 0.5);
            end
        end
    end
end

function plotGridCell(lon, lat, step, color, alpha)
    % Plot a single grid cell
    x = [lon-step/2, lon+step/2, lon+step/2, lon-step/2];
    y = [lat-step/2, lat-step/2, lat+step/2, lat+step/2];
    patch(x, y, color, 'FaceAlpha', alpha, 'EdgeColor', 'none');
end

function plotGroundtrack(longs_deg, lats_deg, longs, colors_array, color, ...
    rp_long, rp_lat, coloured)
    % Plot the groundtrack with optional coloring
    
    % Find discontinuity in longitude
    stop_index = find(diff(longs) > 0, 1);
    
    if ~isempty(stop_index) && stop_index ~= 1
        % Split track
        plotTrackSegment(longs_deg(1:stop_index), lats_deg(1:stop_index), ...
            colors_array(1:stop_index, :), color, coloured);
        plotTrackSegment(longs_deg(stop_index+1:end), lats_deg(stop_index+1:end), ...
            colors_array(stop_index+1:end, :), color, coloured);
    else
        % Continuous track
        plotTrackSegment(longs_deg, lats_deg, colors_array, color, coloured);
    end
    
    % Plot markers
    plot(longs_deg(1), lats_deg(1), 'o', 'MarkerSize', 5, ...
        'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black');
    plot(longs_deg(end), lats_deg(end), 'o', 'MarkerSize', 5, ...
        'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0.929 0.694 0.125]);
    plot(rp_long, rp_lat, 'o', 'MarkerEdgeColor', 'black', ...
        'MarkerFaceColor', 'red', 'MarkerSize', 5);
    
    % Add legend
    addLegend(colors_array, color, coloured);
end

function plotTrackSegment(longs, lats, colors, base_color, coloured)
    % Plot a track segment with optional coloring
    if coloured
        plotColoredSegments(longs, lats, colors);
    else
        plot(longs, lats, '-', 'Color', base_color, 'LineWidth', 2);
    end
end

function plotColoredSegments(longs_deg, lats_deg, colors)
    % Plot line segments with varying colors
    n = length(longs_deg);
    if n < 2, return; end
    
    i = 1;
    while i <= n
        current_color = colors(i, :);
        j = i;
        
        % Find end of current color segment
        while j < n && isequal(colors(j+1, :), current_color)
            j = j + 1;
        end
        
        % Plot segment
        plot(longs_deg(i:j), lats_deg(i:j), '-', ...
            'Color', current_color, 'LineWidth', 2);
        
        % Connect to next segment to avoid gaps
        if j < n
            plot(longs_deg(j:j+1), lats_deg(j:j+1), '-', ...
                'Color', colors(j+1, :), 'LineWidth', 2);
        end
        
        i = j + 1;
    end
end

function addLegend(colors_array, base_color, coloured)
    % Add appropriate legend
    handles = [];
    labels = {};
    
    if coloured
        % Add color legends
        colorMap = {
            base_color, 'Day';
            [0 0 0], 'Night';
            [0.6 0.4 0.2], 'Penumbra';
            [0.4 0.2 0], 'Umbra';
            [1 0 0], 'Invalid Incidence Angle'
        };
        
        for i = 1:size(colorMap, 1)
            if any(all(colors_array == colorMap{i, 1}, 2))
                h = plot(NaN, NaN, '-', 'Color', colorMap{i, 1}, 'LineWidth', 2);
                handles(end+1) = h;
                labels{end+1} = colorMap{i, 2};
            end
        end
    else
        h = plot(NaN, NaN, '-', 'Color', base_color, 'LineWidth', 2);
        handles(end+1) = h;
        labels{end+1} = 'Flyby Groundtrack';
    end
    
    % Add marker legends
    markerInfo = {
        'blue', 'Flyby Entry Point';
        [0.929 0.694 0.125], 'Flyby Exit Point';
        'red', 'Flyby Periapsis'
    };
    
    for i = 1:size(markerInfo, 1)
        h = plot(NaN, NaN, 'o', 'MarkerSize', 5, ...
            'MarkerFaceColor', markerInfo{i, 1}, 'MarkerEdgeColor', 'black');
        handles(end+1) = h;
        labels{end+1} = markerInfo{i, 2};
    end
    
    legend(handles, labels, 'Location', 'best');
end

function incidence_angles = calculateIncidenceAngles(lats, longs, T0)
    % Calculate incidence angles for each point using Enceladus body-fixed frame
    
    % Get positions from Solar System Barycenter (SSB)
    spiceParam.frame = 'J2000';
    spiceParam.abcorr = 'NONE';
    spiceParam.observer = '0';  % SSB
    
    % Get Sun and Enceladus positions in ICRF
    r_sun_icrf = EphSS_car_spice2(10, T0, true, spiceParam);   % Sun position
    r_enc_icrf = EphSS_car_spice2(602, T0, true, spiceParam);  % Enceladus position
    
    % Convert Sun position to Enceladus body-fixed frame
    % This gives us the Sun direction vector in Enceladus body frame
    sun_dir_body = icrf2enceladus(r_sun_icrf, r_enc_icrf, T0);
    
    % Calculate incidence angles for each surface point
    n_points = length(lats);
    incidence_angles = zeros(n_points, 1);
    
    for i = 1:n_points
        % Surface normal at this point in body-fixed frame
        % In body-fixed frame, longitude and latitude directly map to the surface
        normal = [cos(lats(i))*cos(longs(i)); 
                  cos(lats(i))*sin(longs(i)); 
                  sin(lats(i))];
        
        % Calculate incidence angle (0° = overhead, 90° = horizon)
        cos_theta = dot(normal, sun_dir_body);
        incidence_angles(i) = acosd(max(-1, min(1, cos_theta)));
    end
end

% Helper function for inline conditionals
function result = iff(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end