function createAnimatedFlybyPlot(full_state_out_nb, state_LC_resampled, r_enc_history, r_titan_history, full_time_out_nb, pars)
% Creates an interactive, animated plot of the Saturn-centric flyby with
% play/pause, variable speed controls, and direct time input boxes (minutes and hours).
% Now includes both Enceladus and Titan in the animation.
%
% INPUTS:
%   full_state_out_nb      - State history of the N-body trajectory [N x 6]
%   state_LC_resampled     - State history of the Linked Conic trajectory [N x 6]
%   r_enc_history          - Position history of Enceladus [N x 3]
%   r_titan_history        - Position history of Titan [N x 3]
%   full_time_out_nb       - Time vector corresponding to the states [N x 1]
%   pars                   - Parameters struct containing moon and planet data

% --- Animation State Variables ---
isPlaying = false;
playbackSpeed = 1;
playbackDirection = 1; % 1 for forward, -1 for backward
time_in_minutes = full_time_out_nb / 60;
time_in_hours = full_time_out_nb / 3600;
num_steps = length(full_time_out_nb);

% --- Find the initial index for t = -200 minutes ---
[~, initialIndex] = min(abs(time_in_minutes - (-0))); %-200
currentIndex = initialIndex;

% --- Create Figure and Axes ---
fig_anim = figure('Name', 'Animated Saturn-Centric Flyby with Titan and Enceladus', 'Color', 'w', 'Position', [150 150 1000 850]);
ax_anim = axes('Parent', fig_anim, 'Position', [0.1, 0.15, 0.8, 0.75]);
hold(ax_anim, 'on');

% --- Plot Static Elements ---
plot3(ax_anim, 0, 0, 0, 'o', 'MarkerSize', 20, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'DisplayName', 'Saturn');

% --- Add Coordinate axes (length = 2 * Enceladus orbital radius) ---
axis_length = 2 * pars.Moon.OrbRad;
% X-axis arrow (RED)
quiver3(ax_anim, 0, 0, 0, axis_length, 0, 0, 'r', 'LineWidth', 0.5, 'MaxHeadSize', 0.2, 'AutoScale', 'off');
text(ax_anim, axis_length*1.15, 0, 0, 'X', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'r', 'HorizontalAlignment', 'center');
% Y-axis arrow (GREEN)
quiver3(ax_anim, 0, 0, 0, 0, axis_length, 0, 'g', 'LineWidth', 0.5, 'MaxHeadSize', 0.2, 'AutoScale', 'off');
text(ax_anim, 0, axis_length*0.2, 0, 'Y', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'g', 'HorizontalAlignment', 'center');
% Z-axis arrow (BLUE)
quiver3(ax_anim, 0, 0, 0, 0, 0, axis_length, 'b', 'LineWidth', 0.5, 'MaxHeadSize', 0.2, 'AutoScale', 'off');
text(ax_anim, 0, 0, axis_length*1.15, 'Z', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'b', 'HorizontalAlignment', 'center');

plot3(ax_anim, r_enc_history(:,1), r_enc_history(:,2), r_enc_history(:,3), 'k--', 'LineWidth', 1, 'DisplayName', 'Enceladus Orbit');
plot3(ax_anim, r_titan_history(:,1), r_titan_history(:,2), r_titan_history(:,3), 'r--', 'LineWidth', 1, 'DisplayName', 'Titan Orbit');
plot3(ax_anim, full_state_out_nb(:,1), full_state_out_nb(:,2), full_state_out_nb(:,3), ...
      'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.5, 'DisplayName', 'N-Body S/C Trajectory');
plot3(ax_anim, state_LC_resampled(:,1), state_LC_resampled(:,2), state_LC_resampled(:,3), ...
      'Color', [0, 0.4470, 0.7410], 'LineStyle', ':', 'LineWidth', 1.5, 'DisplayName', 'Linked Conic S/C Trajectory');

% --- Initialize Dynamic Plot Handles ---
% Enceladus sphere
radius_enceladus_display = pars.Moon.EquRad * pars.EncPlotSize;
[x_sphere, y_sphere, z_sphere] = sphere(20);
h_enceladus_sphere = surf(ax_anim, [], [], [], 'FaceColor', [0.678, 0.847, 0.902], 'EdgeColor', 'none', 'DisplayName', 'Enceladus', 'FaceAlpha', 0.5);

% Titan sphere (using hardcoded radius of 2574 km)
radius_titan_hardcoded = 2574; % km - Titan's actual radius
radius_titan_display = radius_titan_hardcoded * pars.EncPlotSize/4;
h_titan_sphere = surf(ax_anim, [], [], [], 'FaceColor', [0.8, 0.5, 0.2], 'EdgeColor', 'none', 'DisplayName', 'Titan', 'FaceAlpha', 0.7);

% Spacecraft markers
h_sc_nb = plot3(ax_anim, NaN, NaN, NaN, 'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.8500, 0.3250, 0.0980], 'MarkerEdgeColor', 'k', 'DisplayName', 'N-Body S/C');
h_sc_lc = plot3(ax_anim, NaN, NaN, NaN, 's', 'MarkerSize', 7, 'MarkerFaceColor', [0, 0.4470, 0.7410], 'MarkerEdgeColor', 'k', 'DisplayName', 'Linked Conic S/C');

% --- General Plot Formatting ---
xlabel(ax_anim, 'X [km]');
ylabel(ax_anim, 'Y [km]');
zlabel(ax_anim, 'Z [km]');
axis(ax_anim, 'equal');
grid(ax_anim, 'on');
view(ax_anim, 45, 30);
lighting gouraud;
camlight;
legend(ax_anim, 'show', 'Location', 'northeast');
h_title = title(ax_anim, '');

% --- Create UI Controls ---
h_button_bwd = uicontrol('Parent', fig_anim, 'Style', 'pushbutton', 'String', '<<', ...
                         'Position', [260, 20, 80, 40], 'FontSize', 14, 'Callback', @onBackward);
h_button_play = uicontrol('Parent', fig_anim, 'Style', 'pushbutton', 'String', char(9654), ...
                          'Position', [350, 20, 80, 40], 'FontSize', 14, 'Callback', @onPlayPause);
h_button_fwd = uicontrol('Parent', fig_anim, 'Style', 'pushbutton', 'String', '>>', ...
                         'Position', [440, 20, 80, 40], 'FontSize', 14, 'Callback', @onForward);
h_text_status = uicontrol('Parent', fig_anim, 'Style', 'text', 'String', 'Paused', ...
                          'Position', [290, 70, 200, 25], 'FontSize', 12, 'FontWeight', 'bold');

% Time input controls - Minutes
uicontrol('Parent', fig_anim, 'Style', 'text', 'String', 'Go to Time (min):', ...
          'Position', [550, 60, 120, 20], 'HorizontalAlignment', 'right', 'FontSize', 10);
h_edit_min = uicontrol('Parent', fig_anim, 'Style', 'edit', ...
          'Position', [675, 55, 100, 25], 'FontSize', 10, ...
          'Callback', @onEditTimeMinutes);

% Time input controls - Hours  
uicontrol('Parent', fig_anim, 'Style', 'text', 'String', 'Go to Time (hr):', ...
          'Position', [550, 30, 120, 20], 'HorizontalAlignment', 'right', 'FontSize', 10);
h_edit_hr = uicontrol('Parent', fig_anim, 'Style', 'edit', ...
          'Position', [675, 25, 100, 25], 'FontSize', 10, ...
          'Callback', @onEditTimeHours);

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
            if currentIndex >= num_steps || currentIndex <= 1 % If at an end, restart
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
        if isPlaying, onPlayPause(); end % Pause if playing
        
        val = str2double(get(source, 'String'));
        if isnan(val) || val < time_in_minutes(1) || val > time_in_minutes(end)
            % On invalid input, revert to current time
            set(source, 'String', sprintf('%.2f', time_in_minutes(currentIndex)));
            return;
        end
        
        [~, newIndex] = min(abs(time_in_minutes - val));
        currentIndex = newIndex;
        updateFrame(currentIndex);
    end

    function onEditTimeHours(source, ~)
        if isPlaying, onPlayPause(); end % Pause if playing
        
        val = str2double(get(source, 'String'));
        if isnan(val) || val < time_in_hours(1) || val > time_in_hours(end)
            % On invalid input, revert to current time
            set(source, 'String', sprintf('%.3f', time_in_hours(currentIndex)));
            return;
        end
        
        [~, newIndex] = min(abs(time_in_hours - val));
        currentIndex = newIndex;
        updateFrame(currentIndex);
    end

% --- Helper Functions (Nested) ---
    function updateFrame(k)
        % Updates all visual elements and UI text for a given index k.
        % Update Enceladus sphere
        set(h_enceladus_sphere, 'XData', x_sphere*radius_enceladus_display + r_enc_history(k,1), ...
                              'YData', y_sphere*radius_enceladus_display + r_enc_history(k,2), ...
                              'ZData', z_sphere*radius_enceladus_display + r_enc_history(k,3));
        
        % Update Titan sphere
        set(h_titan_sphere, 'XData', x_sphere*radius_titan_display + r_titan_history(k,1), ...
                          'YData', y_sphere*radius_titan_display + r_titan_history(k,2), ...
                          'ZData', z_sphere*radius_titan_display + r_titan_history(k,3));
        
        % Update spacecraft positions
        set(h_sc_nb, 'XData', full_state_out_nb(k,1), 'YData', full_state_out_nb(k,2), 'ZData', full_state_out_nb(k,3));
        set(h_sc_lc, 'XData', state_LC_resampled(k,1), 'YData', state_LC_resampled(k,2), 'ZData', state_LC_resampled(k,3));
        
        % Update title and time displays
        set(h_title, 'String', sprintf('Saturn-Centric View with Titan & Enceladus | Time: %.2f minutes (%.3f hours)', time_in_minutes(k), time_in_hours(k)));
        set(h_edit_min, 'String', sprintf('%.2f', time_in_minutes(k)));
        set(h_edit_hr, 'String', sprintf('%.3f', time_in_hours(k)));
        drawnow;
    end

    function updateStatusText()
        if ~isPlaying, set(h_text_status, 'String', 'Paused'); return; end
        if playbackDirection == 1, direction_arrow = '>>'; else, direction_arrow = '<<'; end
        set(h_text_status, 'String', sprintf('Speed: %dx %s', playbackSpeed, direction_arrow));
    end

end