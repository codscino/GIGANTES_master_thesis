function liveplot_enceladus_comparison(r_enc_2body, r_enc_nbody, full_time_out_nb, dates, pars)
% Creates an interactive, animated plot comparing 2-body and N-body 
% Enceladus motion with play/pause, variable speed controls, 
% and direct time input via calendar date or MJD2000.
%
% INPUTS:
%   r_enc_2body           - Position history of Enceladus from 2-body approximation [N x 3]
%   r_enc_nbody           - Position history of Enceladus from N-body SPICE [N x 3]
%   full_time_out_nb      - Time vector in seconds from start [N x 1]
%   dates                 - Time vector in MJD2000 [N x 1]
%   pars                  - Parameters struct containing moon and planet data

% --- Animation State Variables ---
isPlaying = false;
playbackSpeed = 1;
playbackDirection = 1; % 1 for forward, -1 for backward
time_in_minutes = full_time_out_nb / 60;
time_in_hours = full_time_out_nb / 3600;
num_steps = length(full_time_out_nb);

% --- Find the initial index for the start time ---
currentIndex = 1;

% --- Create Figure and Axes ---
fig_anim = figure('Name', 'Enceladus Motion: 2-Body vs N-Body Comparison', 'Color', 'w', 'Position', [150 150 1200 900]);
ax_anim = axes('Parent', fig_anim, 'Position', [0.1, 0.15, 0.8, 0.75]);
hold(ax_anim, 'on');

% --- Plot Static Elements ---
plot3(ax_anim, 0, 0, 0, 'o', 'MarkerSize', 25, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'DisplayName', 'Saturn');

% --- Add Coordinate axes (length = 2 * Enceladus orbital radius) ---
axis_length = 2 * pars.Moon.OrbRad;
% X-axis arrow (RED)
quiver3(ax_anim, 0, 0, 0, axis_length, 0, 0, 'r', 'LineWidth', 0.5, 'MaxHeadSize', 0.2, 'AutoScale', 'off');
text(ax_anim, axis_length*1.15, 0, 0, 'X', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'r', 'HorizontalAlignment', 'center');
% Y-axis arrow (GREEN)
quiver3(ax_anim, 0, 0, 0, 0, axis_length, 0, 'g', 'LineWidth', 0.5, 'MaxHeadSize', 0.2, 'AutoScale', 'off');
text(ax_anim, 0, axis_length*1.2, 0, 'Y', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'g', 'HorizontalAlignment', 'center');
% Z-axis arrow (BLUE)
quiver3(ax_anim, 0, 0, 0, 0, 0, axis_length, 'b', 'LineWidth', 0.5, 'MaxHeadSize', 0.2, 'AutoScale', 'off');
text(ax_anim, 0, 0, axis_length*1.15, 'Z', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'b', 'HorizontalAlignment', 'center');

% Plot Enceladus orbital paths
plot3(ax_anim, r_enc_2body(:,1), r_enc_2body(:,2), r_enc_2body(:,3), 'b-', 'LineWidth', 2, 'DisplayName', '2-Body Trajectory');
plot3(ax_anim, r_enc_nbody(:,1), r_enc_nbody(:,2), r_enc_nbody(:,3), 'r-', 'LineWidth', 2, 'DisplayName', 'N-Body (SPICE) Trajectory');

% --- Initialize Dynamic Plot Handles ---
% Enceladus spheres
radius_enceladus_display = pars.Moon.EquRad * pars.EncPlotSize;
[x_sphere, y_sphere, z_sphere] = sphere(20);

% 2-Body Enceladus sphere (Blue)
h_enc_2body = surf(ax_anim, [], [], [], 'FaceColor', 'b', 'EdgeColor', 'none', 'DisplayName', '2-Body Enceladus', 'FaceAlpha', 0.7);

% N-Body Enceladus sphere (Red)
h_enc_nbody = surf(ax_anim, [], [], [], 'FaceColor', 'r', 'EdgeColor', 'none', 'DisplayName', 'N-Body Enceladus', 'FaceAlpha', 0.7);

% Connection line between the two Enceladus positions
h_connection_line = plot3(ax_anim, [NaN, NaN], [NaN, NaN], [NaN, NaN], 'k-', 'LineWidth', 0.5, 'DisplayName', 'Position Difference');

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

% Time input controls - MJD2000
uicontrol('Parent', fig_anim, 'Style', 'text', 'String', 'Go to MJD2000:', ...
          'Position', [550, 60, 120, 20], 'HorizontalAlignment', 'right', 'FontSize', 10);
h_edit_mjd = uicontrol('Parent', fig_anim, 'Style', 'edit', ...
          'Position', [675, 55, 120, 25], 'FontSize', 10, ...
          'Callback', @onEditMJD2000);

% Time input controls - Calendar Date  
uicontrol('Parent', fig_anim, 'Style', 'text', 'String', 'Date (Y-M-D H:M:S):', ...
          'Position', [530, 30, 140, 20], 'HorizontalAlignment', 'right', 'FontSize', 10);
h_edit_date = uicontrol('Parent', fig_anim, 'Style', 'edit', ...
          'Position', [675, 25, 180, 25], 'FontSize', 10, ...
          'Callback', @onEditDateString);

% Statistics display
h_text_stats = uicontrol('Parent', fig_anim, 'Style', 'text', 'String', '', ...
                         'Position', [50, 20, 200, 60], 'FontSize', 9, 'HorizontalAlignment', 'left');

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
            if currentIndex >= num_steps && playbackDirection == 1
                currentIndex = 1; % If at the end, restart from beginning
            elseif currentIndex <= 1 && playbackDirection == -1
                currentIndex = num_steps; % If at the beginning, restart from end
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

    function onEditMJD2000(source, ~)
        if isPlaying, onPlayPause(); end % Pause if playing
        
        val = str2double(get(source, 'String'));
        if isnan(val) || val < dates(1) || val > dates(end)
            % On invalid input, revert to current time
            set(source, 'String', sprintf('%.6f', dates(currentIndex)));
            return;
        end
        
        [~, newIndex] = min(abs(dates - val));
        currentIndex = newIndex;
        updateFrame(currentIndex);
    end

    function onEditDateString(source, ~)
        if isPlaying, onPlayPause(); end % Pause if playing
        
        try
            date_vec = datevec(get(source, 'String'), 'yyyy-mm-dd HH:MM:SS');
            target_mjd = date2mjd2000(date_vec);
            
            if target_mjd < dates(1) || target_mjd > dates(end)
                error('Date is outside the simulation time range.');
            end
            
            [~, newIndex] = min(abs(dates - target_mjd));
            currentIndex = newIndex;
            updateFrame(currentIndex);
        catch ME
            warning('Invalid date format or out of range. Use YYYY-MM-DD HH:MM:SS. Error: %s', ME.message);
            % On invalid input, revert to current time
            current_date_vec = mjd20002date(dates(currentIndex));
            set(source, 'String', datestr(current_date_vec, 'yyyy-mm-dd HH:MM:SS'));
        end
    end

% --- Helper Functions (Nested) ---
    function updateFrame(k)
        % Updates all visual elements and UI text for a given index k.
        
        % Update 2-Body Enceladus sphere (Blue)
        set(h_enc_2body, 'XData', x_sphere*radius_enceladus_display + r_enc_2body(k,1), ...
                         'YData', y_sphere*radius_enceladus_display + r_enc_2body(k,2), ...
                         'ZData', z_sphere*radius_enceladus_display + r_enc_2body(k,3));
        
        % Update N-Body Enceladus sphere (Red)
        set(h_enc_nbody, 'XData', x_sphere*radius_enceladus_display + r_enc_nbody(k,1), ...
                         'YData', y_sphere*radius_enceladus_display + r_enc_nbody(k,2), ...
                         'ZData', z_sphere*radius_enceladus_display + r_enc_nbody(k,3));
        
        % Update connection line
        set(h_connection_line, 'XData', [r_enc_2body(k,1), r_enc_nbody(k,1)], ...
                               'YData', [r_enc_2body(k,2), r_enc_nbody(k,2)], ...
                               'ZData', [r_enc_2body(k,3), r_enc_nbody(k,3)]);
        
        % Calculate current statistics
        curr_diff = norm(r_enc_nbody(k,:) - r_enc_2body(k,:));
        max_diff = max(vecnorm((r_enc_nbody - r_enc_2body)', 2));
        rms_diff = sqrt(mean(vecnorm((r_enc_nbody - r_enc_2body)', 2).^2));
        
        % Get current date vector and format for display
        current_date_vec = mjd20002date(dates(k));
        date_str = datestr(current_date_vec, 'yyyy-mm-dd HH:MM:SS');
        
        % Update title and time displays
        set(h_title, 'String', sprintf('Enceladus: 2-Body vs N-Body | %s | Diff: %.3f km', ...
                                       date_str, curr_diff));
        set(h_edit_mjd, 'String', sprintf('%.6f', dates(k)));
        set(h_edit_date, 'String', date_str);
        
        % Update statistics display
        stats_text = sprintf('Current Diff: %.3f km\nMax Diff: %.3f km\nRMS Diff: %.3f km', curr_diff, max_diff, rms_diff);
        set(h_text_stats, 'String', stats_text);
        
        drawnow;
    end

    function updateStatusText()
        if ~isPlaying, set(h_text_status, 'String', 'Paused'); return; end
        if playbackDirection == 1, direction_arrow = '>>'; else, direction_arrow = '<<'; end
        set(h_text_status, 'String', sprintf('Speed: %dx %s', playbackSpeed, direction_arrow));
    end

end