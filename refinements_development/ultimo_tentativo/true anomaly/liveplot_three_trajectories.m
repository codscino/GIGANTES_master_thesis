function liveplot_three_trajectories(state_LC, state_stm, state_orig, r_enc_history, r_titan_history, full_time_out, pars)
% Creates an interactive, animated plot showing three trajectories:
% 1. Linked Conics trajectory
% 2. STM-corrected N-body trajectory  
% 3. Original N-body trajectory (uncorrected)
%
% INPUTS:
%   state_LC           - Linked Conics trajectory state history [N x 6]
%   state_stm          - STM-corrected N-body trajectory state history [N x 6]
%   state_orig         - Original N-body trajectory state history [N x 6]
%   r_enc_history      - Position history of Enceladus [N x 3]
%   r_titan_history    - Position history of Titan [N x 3]
%   full_time_out      - Time vector in seconds [N x 1]
%   pars               - Parameters struct containing moon and planet data

% --- Animation State Variables ---
isPlaying = false;
playbackSpeed = 1;
playbackDirection = 1; % 1 for forward, -1 for backward
time_in_minutes = full_time_out / 60;
time_in_hours = full_time_out / 3600;
num_steps = length(full_time_out);

% --- Planet Size Control ---
EncPlotSize = 1; % Default to normal size

% --- STM Following Control ---
followSTM = false; % Default to not following STM

% --- Pre-calculate True Anomalies for all time steps ---
true_anomalies_deg = zeros(num_steps, 1);
eccentricities = zeros(num_steps, 1);
valid_indices = false(num_steps, 1);


for i = 1:num_steps
    try
        kep_stm = car2kep(state_stm(i,:), pars.Planet.mu);
        true_anomalies_deg(i) = rad2deg(kep_stm(6));
        eccentricities(i) = kep_stm(2);
        valid_indices(i) = true;
    catch
        % Keep NaN values for invalid orbital elements
        true_anomalies_deg(i) = NaN;
        eccentricities(i) = NaN;
        valid_indices(i) = false;
    end
end

% --- Find the initial index for t = -2 minutes ---
[~, initialIndex] = min(abs(time_in_minutes + 2)); % Changed from 0 to -2 minutes
currentIndex = initialIndex;

% --- Create Figure and Axes ---
fig_title = 'Three-Trajectory Comparison: Linked Conics vs STM-Corrected N-Body vs Original N-Body';
fig_anim = figure('Name', fig_title, 'Color', 'w', 'Position', [100 100 1200 900]);
ax_anim = axes('Parent', fig_anim, 'Position', [0.1, 0.15, 0.8, 0.75]);
hold(ax_anim, 'on');

% --- Plot Static Elements ---
plot3(ax_anim, 0, 0, 0, 'o', 'MarkerSize', 20, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'DisplayName', 'Saturn');


% --- Add Coordinate axes ---
axis_length = 2 * pars.Moon.OrbRad;
% X-axis arrow (RED)
quiver3(ax_anim, 0, 0, 0, axis_length, 0, 0, 'r', 'LineWidth', 0.5, 'MaxHeadSize', 0.2, 'AutoScale', 'off', 'HandleVisibility', 'off');
text(ax_anim, axis_length*1.15, 0, 0, 'X', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'r', 'HorizontalAlignment', 'center');
% Y-axis arrow (GREEN)
quiver3(ax_anim, 0, 0, 0, 0, axis_length, 0, 'g', 'LineWidth', 0.5, 'MaxHeadSize', 0.2, 'AutoScale', 'off', 'HandleVisibility', 'off');
text(ax_anim, 0, axis_length*0.2, 0, 'Y', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'g', 'HorizontalAlignment', 'center');
% Z-axis arrow (BLUE)
quiver3(ax_anim, 0, 0, 0, 0, 0, axis_length, 'b', 'LineWidth', 0.5, 'MaxHeadSize', 0.2, 'AutoScale', 'off', 'HandleVisibility', 'off');
text(ax_anim, 0, 0, axis_length*1.15, 'Z', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'b', 'HorizontalAlignment', 'center');
% --- Plot Orbital Paths and Trajectories ---
plot3(ax_anim, r_enc_history(:,1), r_enc_history(:,2), r_enc_history(:,3), 'k--', 'LineWidth', 1, 'DisplayName', 'Enceladus Orbit');
plot3(ax_anim, r_titan_history(:,1), r_titan_history(:,2), r_titan_history(:,3), 'r--', 'LineWidth', 1, 'DisplayName', 'Titan Orbit');

% Plot the three spacecraft trajectories
plot3(ax_anim, state_LC(:,1), state_LC(:,2), state_LC(:,3), ...
      'Color', [0, 0.4470, 0.7410 0.2], 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Linked Conics S/C');
plot3(ax_anim, state_stm(:,1), state_stm(:,2), state_stm(:,3), ...
      'Color', [0.8500, 0.3250, 0.0980 0.8], 'LineWidth', 1, 'DisplayName', 'STM-Corrected N-Body S/C');
plot3(ax_anim, state_orig(:,1), state_orig(:,2), state_orig(:,3), ...
      'Color', [0, 1, 0 0.2], 'LineStyle', ':', 'LineWidth', 3, 'DisplayName', 'Original N-Body S/C');

% --- Initialize Dynamic Plot Handles ---
% Enceladus sphere
[x_sphere, y_sphere, z_sphere] = sphere(20);
h_enceladus_sphere = surf(ax_anim, [], [], [], 'FaceColor', [0.678, 0.847, 0.902], 'EdgeColor', 'none', 'DisplayName', 'Enceladus', 'FaceAlpha', 0.5);

% Titan sphere
radius_titan_hardcoded = 2574; % km - Titan's actual radius
h_titan_sphere = surf(ax_anim, [], [], [], 'FaceColor', [0.8, 0.5, 0.2], 'EdgeColor', 'none', 'DisplayName', 'Titan', 'FaceAlpha', 0.7);

% Spacecraft markers for all three trajectories
h_sc_lc = plot3(ax_anim, NaN, NaN, NaN, 's', 'MarkerSize', 8, 'MarkerFaceColor', [0, 0.4470, 0.7410], 'MarkerEdgeColor', 'k', 'DisplayName', 'LC S/C Position');
h_sc_stm = plot3(ax_anim, NaN, NaN, NaN, 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.8500, 0.3250, 0.0980], 'MarkerEdgeColor', 'k', 'DisplayName', 'STM S/C Position');
h_sc_orig = plot3(ax_anim, NaN, NaN, NaN, '^', 'MarkerSize', 8, 'MarkerFaceColor', [0.9290, 0.6940, 0.1250], 'MarkerEdgeColor', 'k', 'DisplayName', 'Orig S/C Position');

% --- General Plot Formatting ---
xlabel(ax_anim, 'X [km]');
ylabel(ax_anim, 'Y [km]');
zlabel(ax_anim, 'Z [km]');
axis(ax_anim, 'equal');
grid(ax_anim, 'on');
view(ax_anim, 45, 30);

% --- Set Initial Zoom on Enceladus (DEFAULT) ---
% Get Enceladus position at initial time (t=-2 minutes)
enc_pos_initial = r_enc_history(initialIndex, :);

% Define zoom area as 3 times Enceladus radius (new default)
zoom_area = 3 * pars.Moon.EquRad;

% Set axis limits centered on Enceladus with square area
xlim(ax_anim, [enc_pos_initial(1) - zoom_area, enc_pos_initial(1) + zoom_area]);
ylim(ax_anim, [enc_pos_initial(2) - zoom_area, enc_pos_initial(2) + zoom_area]);
zlim(ax_anim, [enc_pos_initial(3) - zoom_area, enc_pos_initial(3) + zoom_area]);

% Store the full plot limits for later use (calculate from all trajectory data)
all_x_data = [state_LC(:,1); state_stm(:,1); state_orig(:,1); r_enc_history(:,1); r_titan_history(:,1)];
all_y_data = [state_LC(:,2); state_stm(:,2); state_orig(:,2); r_enc_history(:,2); r_titan_history(:,2)];
all_z_data = [state_LC(:,3); state_stm(:,3); state_orig(:,3); r_enc_history(:,3); r_titan_history(:,3)];

full_xlim = [min(all_x_data), max(all_x_data)];
full_ylim = [min(all_y_data), max(all_y_data)];
full_zlim = [min(all_z_data), max(all_z_data)];

% Add some padding to full limits
padding = 0.1 * max([max(all_x_data)-min(all_x_data), max(all_y_data)-min(all_y_data), max(all_z_data)-min(all_z_data)]);
full_xlim = full_xlim + [-padding, padding];
full_ylim = full_ylim + [-padding, padding];
full_zlim = full_zlim + [-padding, padding];

lighting gouraud;
camlight;
legend(ax_anim, 'show', 'Location', 'northeast');
h_title = title(ax_anim, '');

% Disable clipping to prevent objects from disappearing when zooming
set(ax_anim, 'Clipping', 'off');

% --- Create UI Controls ---
h_button_bwd = uicontrol('Parent', fig_anim, 'Style', 'pushbutton', 'String', '<<', ...
                         'Position', [260, 20, 80, 40], 'FontSize', 14, 'Callback', @onBackward);
h_button_play = uicontrol('Parent', fig_anim, 'Style', 'pushbutton', 'String', char(9654), ...
                          'Position', [350, 20, 80, 40], 'FontSize', 14, 'Callback', @onPlayPause);
h_button_fwd = uicontrol('Parent', fig_anim, 'Style', 'pushbutton', 'String', '>>', ...
                         'Position', [440, 20, 80, 40], 'FontSize', 14, 'Callback', @onForward);
h_text_status = uicontrol('Parent', fig_anim, 'Style', 'text', 'String', 'Paused', ...
                          'Position', [290, 70, 200, 25], 'FontSize', 12, 'FontWeight', 'bold');

% --- NEW: Three Zoom Control Buttons ---
h_button_zoom_enc = uicontrol('Parent', fig_anim, 'Style', 'pushbutton', 'String', 'Enceladus Zoom', ...
                             'Position', [750, 100, 100, 30], 'FontSize', 9, 'Callback', @onZoomEnceladus, ...
                             'BackgroundColor', [0.9, 1, 0.9]); % Light green to indicate default

h_button_zoom_flyby = uicontrol('Parent', fig_anim, 'Style', 'pushbutton', 'String', 'Flyby Zoom', ...
                               'Position', [860, 100, 100, 30], 'FontSize', 9, 'Callback', @onZoomFlyby);

h_button_zoom_full = uicontrol('Parent', fig_anim, 'Style', 'pushbutton', 'String', 'Full View', ...
                              'Position', [970, 100, 80, 30], 'FontSize', 9, 'Callback', @onZoomFull);

% Large planets toggle
h_checkbox_large = uicontrol('Parent', fig_anim, 'Style', 'checkbox', 'String', 'Large Planets', ...
                            'Position', [1060, 100, 100, 30], 'FontSize', 10, 'Value', 0, ...
                            'Callback', @onToggleLargePlanets);

% Follow STM button
h_button_follow_stm = uicontrol('Parent', fig_anim, 'Style', 'pushbutton', 'String', 'Follow STM', ...
                               'Position', [1060, 130, 80, 30], 'FontSize', 9, 'Callback', @onToggleFollowSTM);

% Time input controls - Minutes
uicontrol('Parent', fig_anim, 'Style', 'text', 'String', 'Go to Time (min):', ...
          'Position', [550, 90, 120, 20], 'HorizontalAlignment', 'right', 'FontSize', 10);
h_edit_min = uicontrol('Parent', fig_anim, 'Style', 'edit', ...
          'Position', [675, 85, 100, 25], 'FontSize', 10, ...
          'Callback', @onEditTimeMinutes);

% Time input controls - Hours  
uicontrol('Parent', fig_anim, 'Style', 'text', 'String', 'Go to Time (hr):', ...
          'Position', [550, 60, 120, 20], 'HorizontalAlignment', 'right', 'FontSize', 10);
h_edit_hr = uicontrol('Parent', fig_anim, 'Style', 'edit', ...
          'Position', [675, 55, 100, 25], 'FontSize', 10, ...
          'Callback', @onEditTimeHours);

% NEW: True Anomaly input control
uicontrol('Parent', fig_anim, 'Style', 'text', 'String', 'Go to True Anom (°):', ...
          'Position', [550, 30, 120, 20], 'HorizontalAlignment', 'right', 'FontSize', 10);
h_edit_trueanom = uicontrol('Parent', fig_anim, 'Style', 'edit', ...
          'Position', [675, 25, 100, 25], 'FontSize', 10, ...
          'Callback', @onEditTrueAnomaly, 'BackgroundColor', [1, 1, 0.9]);

% --- Distance Comparison Display ---
uicontrol('Parent', fig_anim, 'Style', 'text', 'String', 'Trajectory Separations:', ...
          'Position', [20, 90, 200, 20], 'HorizontalAlignment', 'left', 'FontSize', 11, 'FontWeight', 'bold');

uicontrol('Parent', fig_anim, 'Style', 'text', 'String', 'LC vs STM:', ...
          'Position', [20, 70, 80, 15], 'HorizontalAlignment', 'left', 'FontSize', 9);
h_text_dist_lc_stm = uicontrol('Parent', fig_anim, 'Style', 'text', 'String', '0.000 km', ...
                               'Position', [100, 70, 80, 15], 'FontSize', 9, ...
                               'BackgroundColor', [0.95, 0.95, 0.95], 'HorizontalAlignment', 'center');

uicontrol('Parent', fig_anim, 'Style', 'text', 'String', 'LC vs Orig:', ...
          'Position', [20, 50, 80, 15], 'HorizontalAlignment', 'left', 'FontSize', 9);
h_text_dist_lc_orig = uicontrol('Parent', fig_anim, 'Style', 'text', 'String', '0.000 km', ...
                                'Position', [100, 50, 80, 15], 'FontSize', 9, ...
                                'BackgroundColor', [0.95, 0.95, 0.95], 'HorizontalAlignment', 'center');

uicontrol('Parent', fig_anim, 'Style', 'text', 'String', 'STM vs Orig:', ...
          'Position', [20, 30, 80, 15], 'HorizontalAlignment', 'left', 'FontSize', 9);
h_text_dist_stm_orig = uicontrol('Parent', fig_anim, 'Style', 'text', 'String', '0.000 km', ...
                                 'Position', [100, 30, 80, 15], 'FontSize', 9, ...
                                 'BackgroundColor', [0.95, 0.95, 0.95], 'HorizontalAlignment', 'center');

% --- Initialize Plot and UI to Default State ---
updateFrame(currentIndex);

% --- Main Animation Loop (Nested Function) ---
    function runAnimationLoop()
        while isPlaying
            if ~isvalid(fig_anim), break; end
            
            currentIndex = currentIndex + (playbackSpeed * playbackDirection);
            
            if currentIndex >= num_steps || currentIndex <= 1
                currentIndex = max(1, min(num_steps, currentIndex));
                isPlaying = false;
            end
            
            updateFrame(round(currentIndex));
            
            if ~isPlaying
                set(h_button_play, 'String', char(9654));
                updateStatusText();
            end
            
            pause(0.01);
        end
    end

% --- UI Callback Functions (Nested) ---
    function onPlayPause(~,~)
        isPlaying = ~isPlaying;
        if isPlaying
            set(h_button_play, 'String', char(9616));
            if currentIndex >= num_steps || currentIndex <= 1
                currentIndex = 1;
                playbackDirection = 1;
                playbackSpeed = 1;
            end
        else
            set(h_button_play, 'String', char(9654));
        end
        updateStatusText();
        
        if isPlaying, runAnimationLoop(); end
    end

    function onForward(~,~)
        if playbackDirection == -1, playbackSpeed = 1; else, playbackSpeed = min(1024, playbackSpeed * 2); end
        playbackDirection = 1;
        updateStatusText();
        if ~isPlaying, onPlayPause(); end
    end

    function onBackward(~,~)
        if playbackDirection == 1, playbackSpeed = 1; else, playbackSpeed = min(1024, playbackSpeed * 2); end
        playbackDirection = -1;
        updateStatusText();
        if ~isPlaying, onPlayPause(); end
    end

    function onEditTimeMinutes(source, ~)
        if isPlaying, onPlayPause(); end
        
        val = str2double(get(source, 'String'));
        if isnan(val) || val < time_in_minutes(1) || val > time_in_minutes(end)
            set(source, 'String', sprintf('%.2f', time_in_minutes(currentIndex)));
            return;
        end
        
        [~, newIndex] = min(abs(time_in_minutes - val));
        currentIndex = newIndex;
        updateFrame(currentIndex);
    end

    function onEditTimeHours(source, ~)
        if isPlaying, onPlayPause(); end
        
        val = str2double(get(source, 'String'));
        if isnan(val) || val < time_in_hours(1) || val > time_in_hours(end)
            set(source, 'String', sprintf('%.3f', time_in_hours(currentIndex)));
            return;
        end
        
        [~, newIndex] = min(abs(time_in_hours - val));
        currentIndex = newIndex;
        updateFrame(currentIndex);
    end

    % NEW: True Anomaly input callback
    function onEditTrueAnomaly(source, ~)
        if isPlaying, onPlayPause(); end
        
        val = str2double(get(source, 'String'));
        if isnan(val)
            % Reset to current value if invalid input
            if valid_indices(currentIndex)
                set(source, 'String', sprintf('%.3f', true_anomalies_deg(currentIndex)));
            else
                set(source, 'String', 'N/A');
            end
            return;
        end
        
        % Normalize input to 0-360 range
        val = mod(val, 360);
        
        % Find closest valid true anomaly
        valid_ta = true_anomalies_deg(valid_indices);
        valid_idx_map = find(valid_indices);
        
        if isempty(valid_ta)
            fprintf('No valid orbital elements found in trajectory!\n');
            return;
        end
        
        % Handle wrap-around by considering both val and val±360
        candidates = [valid_ta; valid_ta + 360; valid_ta - 360];
        candidate_indices = [valid_idx_map; valid_idx_map; valid_idx_map];
        
        [~, min_idx] = min(abs(candidates - val));
        newIndex = candidate_indices(min_idx);
        
        % Ensure index is within bounds
        newIndex = max(1, min(num_steps, newIndex));
        
        currentIndex = newIndex;
        updateFrame(currentIndex);
  
    end

    % --- NEW: Three Zoom Functions ---
    function onZoomEnceladus(~,~)
        % Disable STM following
        followSTM = false;
        set(h_button_follow_stm, 'String', 'Follow STM', 'BackgroundColor', [0.94, 0.94, 0.94]);
        
        % Get current Enceladus position
        enc_pos_current = r_enc_history(currentIndex, :);
        zoom_area = 3 * pars.Moon.EquRad; % Close-up view
        
        xlim(ax_anim, [enc_pos_current(1) - zoom_area, enc_pos_current(1) + zoom_area]);
        ylim(ax_anim, [enc_pos_current(2) - zoom_area, enc_pos_current(2) + zoom_area]);
        zlim(ax_anim, [enc_pos_current(3) - zoom_area, enc_pos_current(3) + zoom_area]);
        
        % Visual feedback - highlight current button
        set(h_button_zoom_enc, 'BackgroundColor', [0.9, 1, 0.9]);
        set(h_button_zoom_flyby, 'BackgroundColor', [0.94, 0.94, 0.94]);
        set(h_button_zoom_full, 'BackgroundColor', [0.94, 0.94, 0.94]);
    end

    function onZoomFlyby(~,~)
        % Disable STM following
        followSTM = false;
        set(h_button_follow_stm, 'String', 'Follow STM', 'BackgroundColor', [0.94, 0.94, 0.94]);
        
        % Get current Enceladus position
        enc_pos_current = r_enc_history(currentIndex, :);
        zoom_area = 30 * pars.Moon.EquRad; % Flyby context view
        
        xlim(ax_anim, [enc_pos_current(1) - zoom_area, enc_pos_current(1) + zoom_area]);
        ylim(ax_anim, [enc_pos_current(2) - zoom_area, enc_pos_current(2) + zoom_area]);
        zlim(ax_anim, [enc_pos_current(3) - zoom_area, enc_pos_current(3) + zoom_area]);
        
        % Visual feedback - highlight current button
        set(h_button_zoom_enc, 'BackgroundColor', [0.94, 0.94, 0.94]);
        set(h_button_zoom_flyby, 'BackgroundColor', [0.9, 1, 0.9]);
        set(h_button_zoom_full, 'BackgroundColor', [0.94, 0.94, 0.94]);
    end

    function onZoomFull(~,~)
        % Disable STM following
        followSTM = false;
        set(h_button_follow_stm, 'String', 'Follow STM', 'BackgroundColor', [0.94, 0.94, 0.94]);
        
        xlim(ax_anim, full_xlim);
        ylim(ax_anim, full_ylim);
        zlim(ax_anim, full_zlim);
        
        % Visual feedback - highlight current button
        set(h_button_zoom_enc, 'BackgroundColor', [0.94, 0.94, 0.94]);
        set(h_button_zoom_flyby, 'BackgroundColor', [0.94, 0.94, 0.94]);
        set(h_button_zoom_full, 'BackgroundColor', [0.9, 1, 0.9]);
    end

    function onToggleLargePlanets(source, ~)
        % Update EncPlotSize based on checkbox state
        if get(source, 'Value')
            EncPlotSize = 70; % Large planets
        else
            EncPlotSize = 1;  % Normal size
        end
        
        % Refresh current frame to show updated sizes
        updateFrame(currentIndex);
    end

    function onToggleFollowSTM(~,~)
        followSTM = ~followSTM;
        if followSTM
            set(h_button_follow_stm, 'String', 'Following STM', 'BackgroundColor', [0.9, 1, 0.9]);
            % Reset other zoom button highlights
            set(h_button_zoom_enc, 'BackgroundColor', [0.94, 0.94, 0.94]);
            set(h_button_zoom_flyby, 'BackgroundColor', [0.94, 0.94, 0.94]);
            set(h_button_zoom_full, 'BackgroundColor', [0.94, 0.94, 0.94]);
            % Immediately apply STM following
            updateSTMFollowView(currentIndex);
        else
            set(h_button_follow_stm, 'String', 'Follow STM', 'BackgroundColor', [0.94, 0.94, 0.94]);
        end
    end

% --- Helper Functions (Nested) ---
    function updateSTMFollowView(k)
        % Update view to follow STM spacecraft position
        stm_pos_current = state_stm(k,1:3);
        zoom_area = 3 * pars.Moon.EquRad; % Same as Enceladus zoom
        
        xlim(ax_anim, [stm_pos_current(1) - zoom_area, stm_pos_current(1) + zoom_area]);
        ylim(ax_anim, [stm_pos_current(2) - zoom_area, stm_pos_current(2) + zoom_area]);
        zlim(ax_anim, [stm_pos_current(3) - zoom_area, stm_pos_current(3) + zoom_area]);
    end

    function updateFrame(k)
        % Recalculate sphere sizes based on current EncPlotSize
        current_radius_enceladus = pars.Moon.EquRad * EncPlotSize;
        current_radius_titan = radius_titan_hardcoded * EncPlotSize/5;
        
        % Update Enceladus sphere
        set(h_enceladus_sphere, 'XData', x_sphere*current_radius_enceladus + r_enc_history(k,1), ...
                              'YData', y_sphere*current_radius_enceladus + r_enc_history(k,2), ...
                              'ZData', z_sphere*current_radius_enceladus + r_enc_history(k,3));
        
        % Update Titan sphere
        set(h_titan_sphere, 'XData', x_sphere*current_radius_titan + r_titan_history(k,1), ...
                          'YData', y_sphere*current_radius_titan + r_titan_history(k,2), ...
                          'ZData', z_sphere*current_radius_titan + r_titan_history(k,3));
        
        % Update all three spacecraft positions
        set(h_sc_lc, 'XData', state_LC(k,1), 'YData', state_LC(k,2), 'ZData', state_LC(k,3));
        set(h_sc_stm, 'XData', state_stm(k,1), 'YData', state_stm(k,2), 'ZData', state_stm(k,3));
        set(h_sc_orig, 'XData', state_orig(k,1), 'YData', state_orig(k,2), 'ZData', state_orig(k,3));
        
        % Calculate trajectory separations
        pos_lc = state_LC(k,1:3);
        pos_stm = state_stm(k,1:3);
        pos_orig = state_orig(k,1:3);
        
        dist_lc_stm = norm(pos_lc - pos_stm);
        dist_lc_orig = norm(pos_lc - pos_orig);
        dist_stm_orig = norm(pos_stm - pos_orig);
        
        set(h_text_dist_lc_stm, 'String', sprintf('%.3f km', dist_lc_stm));
        set(h_text_dist_lc_orig, 'String', sprintf('%.3f km', dist_lc_orig));
        set(h_text_dist_stm_orig, 'String', sprintf('%.3f km', dist_stm_orig));
        
        % Update view if following STM
        if followSTM
            updateSTMFollowView(k);
        end
        
        % Update title and time displays
        set(h_title, 'String', sprintf('Three-Trajectory Comparison | Time: %.2f minutes (%.3f hours)', time_in_minutes(k), time_in_hours(k)));
        set(h_edit_min, 'String', sprintf('%.2f', time_in_minutes(k)));
        set(h_edit_hr, 'String', sprintf('%.3f', time_in_hours(k)));
        
        % Update true anomaly input display
        if valid_indices(k)
            set(h_edit_trueanom, 'String', sprintf('%.3f', true_anomalies_deg(k)));
        else
            set(h_edit_trueanom, 'String', 'N/A');
        end
        
        drawnow;
    end

    function updateStatusText()
        if ~isPlaying, set(h_text_status, 'String', 'Paused'); return; end
        if playbackDirection == 1, direction_arrow = '>>'; else, direction_arrow = '<<'; end
        set(h_text_status, 'String', sprintf('Speed: %dx %s', playbackSpeed, direction_arrow));
    end

end