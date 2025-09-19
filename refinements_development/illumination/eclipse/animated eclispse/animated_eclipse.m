function animated_eclipse(T_start, T_end, time_step_seconds, grid_step_deg, use_parallel)
%
% Animates the lighting conditions on Enceladus during Eclipse,
% visualizing areas under umbra and penumbra using a scatter plot approach.
% The animation pre-computes all frames and provides GUI controls for
% playback and scrubbing.
%
%   INPUTS:
%   T_start             - Start time for the animation in MJD2000 days
%   T_end               - End time for the animation in MJD2000 days
%   time_step_seconds   - Time step between animation frames in seconds
%   grid_step_deg       - illumination grid resolution in degrees
%   use_parallel        - Flag to enable or disable parallel pre-computation (boolean).
%                         If not provided, defaults to 'false'.


% Handle default input for use_parallel flag
if nargin < 5
    use_parallel = false;
end



%% 1. Initialization and Setup
disp('Initializing parameters and loading SPICE kernels...');
kernels = {'sat441.bsp', 'naif0012.tls'};
loadSpiceKernels(kernels);
R_enc = 252.1;
spiceParam.frame = 'J2000';
spiceParam.abcorr = 'NONE';
spiceParam.observer = '602';

%% 2. Define Time Range and Surface Grid from Inputs
seconds_in_a_day = 86400;
T_step = time_step_seconds / seconds_in_a_day; % Convert seconds to days
time_vector_mjd = T_start:T_step:T_end;
num_time_steps = length(time_vector_mjd);
fprintf('Time range defined for %d steps.\n', num_time_steps);

grid_step_rad = deg2rad(grid_step_deg);
lat_vector = (-pi/2:grid_step_rad:pi/2)';
lon_vector1 = (0:grid_step_rad:pi/2)';
lon_vector2 = (3/2*pi:grid_step_rad:2*pi)';
lon_vector = [lon_vector1;lon_vector2];
[lat_grid, lon_grid] = meshgrid(lat_vector, lon_vector);
lat_array = lat_grid(:);
lon_array = lon_grid(:);
num_points = numel(lat_array);

%% 3. Setup Figure and Axes for Pre-computation
disp('Setting up the main figure...');
fig = figure('Name', 'Enceladus Lighting Animation (Scatter)', 'Position', [100, 100, 900, 750], 'Visible', 'off');
ax = axes('Parent', fig, 'Position', [0.1, 0.25, 0.8, 0.65]);
plotTextureLatLong_revised('Enceladus', ax);
hold(ax, 'on');

% --- FIX #1: Lock the axis limits to the map boundaries ---
ylim(ax, [-90 90]);
xlim(ax, [0 360]);

% --- FIX #3: Enable clipping to hide marker parts outside boundaries ---
set(ax, 'Clipping', 'on');
set(ax, 'ClippingStyle', 'rectangle');
% -----------------------------------------------------------

title(ax, 'Loading... Please wait for pre-computation to finish.');

% --- Start Parallel Pool and Initialize Workers ---
if use_parallel
  loadParallelKernels(kernels)
end

%% 4. Pre-compute and Pre-render All Plot Objects
disp('Pre-calculating and rendering all frames using scatter approach...');
plot_handles_all_steps = cell(1, num_time_steps);
progress_bar = waitbar(0, 'Pre-rendering all animation frames...');

for i = 1:num_time_steps
    T = time_vector_mjd(i);
    T_array = T * ones(num_points, 1);
    dataTable = table(lat_array, lon_array, T_array, 'VariableNames', {'lat', 'lon', 'T'});

    lighting_status = calculate_lighting_conditions(dataTable, R_enc, spiceParam, use_parallel);

    dataTable.Lighting = lighting_status;

    frame_handles = createBoundaryPatchesScatter(dataTable, ax, grid_step_deg);

    plot_handles_all_steps{i} = frame_handles;
    waitbar(i / num_time_steps, progress_bar);
end
close(progress_bar);
hold(ax, 'off');
disp('All frames have been rendered.');

% =========================================================================
% ======================= NEW CODE FOR LEGEND STARTS HERE =======================
% =========================================================================
%% 5. Create Legend
disp('Creating plot legend...');
hold(ax, 'on');

% Create dummy scatter plots with NaN coordinates to generate handles for the legend
% without plotting any visible data points on the main axes.
% The properties (color, marker, alpha) must match the actual data plots.
h_umbra_legend = scatter(ax, NaN, NaN, 100, 'black', 'filled', 's', ...
    'DisplayName', 'Umbra', 'MarkerFaceAlpha', 0.4);
h_penumbra_legend = scatter(ax, NaN, NaN, 100, [0.7 0.7 0.7], 'filled', 's', ...
    'DisplayName', 'Penumbra', 'MarkerFaceAlpha', 0.4);

% Create the legend using the handles of the dummy plots
lgd = legend([h_umbra_legend, h_penumbra_legend]);
lgd.Location = 'northwest'; % Place in the top-left corner
lgd.Color = [0.9 0.9 0.9];  % Set a light gray background color for readability
lgd.TextColor = 'black';   % Ensure text is readable
lgd.FontSize = 10;
title(lgd, 'Eclipse Type'); % Add a title to the legend

hold(ax, 'off');
% =========================================================================
% ======================== NEW CODE FOR LEGEND ENDS HERE ========================
% =========================================================================

%% 6. Create GUI Controls (Unchanged)
time_label = uicontrol('Parent', fig, 'Style', 'text',...
    'String', datestr(mjd20002datenum(time_vector_mjd(1))),...
    'Position', [350, 90, 200, 25],...
    'FontSize', 14);

slider = uicontrol('Parent', fig, 'Style', 'slider',...
    'Min', 1, 'Max', num_time_steps, 'Value', 1,...
    'Position', [100, 60, 700, 20],...
    'SliderStep', [1/(num_time_steps-1), 10/(num_time_steps-1)], ...
    'Callback', @slider_callback);

play_button = uicontrol('Parent', fig, 'Style', 'pushbutton',...
    'String', 'Play',...
    'Position', [400, 10, 100, 40],...
    'FontSize', 12,...
    'Callback', @play_pause_callback);

play_button.UserData = struct('is_playing', false);

%% 7. Final Setup and Callback Definitions
current_index = 1;
previous_index = -1;
updatePlot(1);
set(fig, 'Visible', 'on');
title(ax, 'Enceladus Lighting Conditions');
disp('Animation is ready.');

% --- Callback Functions ---
function slider_callback(source, ~)
    new_index = round(source.Value);
    updatePlot(new_index);
end

function play_pause_callback(source, ~)
    state = source.UserData;
    state.is_playing = ~state.is_playing;
    source.UserData = state;
    
    if state.is_playing
        source.String = 'Pause';
        runAnimation();
    else
        source.String = 'Play';
    end
end

function runAnimation()
    while play_button.UserData.is_playing && isvalid(play_button)
        current_index = current_index + 1;
        if current_index > num_time_steps
            current_index = 1;
        end
        updatePlot(current_index);
        slider.Value = current_index;
        pause(0.05);
        drawnow;
    end
    if isvalid(play_button)
        play_button.UserData.is_playing = false;
        play_button.String = 'Play';
    end
end

function updatePlot(new_index)
    if new_index == previous_index
        return;
    end
    
    if previous_index > 0 && previous_index <= num_time_steps
        handles_to_hide = plot_handles_all_steps{previous_index};
        if isfield(handles_to_hide, 'penumbra') && ~isempty(handles_to_hide.penumbra) && isvalid(handles_to_hide.penumbra)
            set(handles_to_hide.penumbra, 'Visible', 'off');
        end
        if isfield(handles_to_hide, 'umbra') && ~isempty(handles_to_hide.umbra) && isvalid(handles_to_hide.umbra)
            set(handles_to_hide.umbra, 'Visible', 'off');
        end
    end
    
    handles_to_show = plot_handles_all_steps{new_index};
    if isfield(handles_to_show, 'penumbra') && ~isempty(handles_to_show.penumbra) && isvalid(handles_to_show.penumbra)
        set(handles_to_show.penumbra, 'Visible', 'on');
    end
    if isfield(handles_to_show, 'umbra') && ~isempty(handles_to_show.umbra) && isvalid(handles_to_show.umbra)
        set(handles_to_show.umbra, 'Visible', 'on');
    end
    
    current_datenum = mjd20002datenum(time_vector_mjd(new_index));
    set(time_label, 'String', datestr(current_datenum, 'yyyy-mm-dd HH:MM:SS'));
    
    current_index = new_index;
    previous_index = new_index;
end
%% Helper functions
function dn = mjd20002datenum(mjd)
    dn = mjd + datenum('2000-01-01 12:00:00');
end

function frame_handles = createBoundaryPatchesScatter(dataTable, ax, step_deg)
    frame_handles = struct('penumbra', [], 'umbra', []);
    
    lat_deg = rad2deg(dataTable.lat);
    lon_deg = rad2deg(dataTable.lon);
    
    marker_size = (step_deg*3)^2; 
    
    valid_indices = (lon_deg <= 90) | (lon_deg >= 270);
    
    penumbra_indices = find(strcmp(dataTable.Lighting, 'Penumbra') & valid_indices);
    if ~isempty(penumbra_indices)
        h = scatter(ax, lon_deg(penumbra_indices), lat_deg(penumbra_indices), ...
                   marker_size, [0.7 0.7 0.7], 'filled', 's', ...
                   'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none', ...
                   'Visible', 'off');
        set(h, 'Clipping', 'on');
        frame_handles.penumbra = h;
    end
    
    umbra_indices = find(strcmp(dataTable.Lighting, 'Umbra') & valid_indices);
    if ~isempty(umbra_indices)
        h = scatter(ax, lon_deg(umbra_indices), lat_deg(umbra_indices), ...
                   marker_size, 'black', 'filled', 's', ...
                   'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none', ...
                   'Visible', 'off');
        set(h, 'Clipping', 'on');
        frame_handles.umbra = h;
    end
end
end