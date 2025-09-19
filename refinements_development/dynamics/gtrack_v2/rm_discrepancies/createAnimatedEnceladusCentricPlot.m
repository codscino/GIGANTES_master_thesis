function createAnimatedEnceladusCentricPlot(r_sc_in_flyby_frame_nb, r_sc_in_flyby_frame_LC, r_saturn_in_flyby_frame_nb, full_time_out_nb, pars)
% Creates an interactive, animated plot of the Enceladus-centric flyby
% with play/pause, variable speed controls, and direct time input boxes (minutes and hours).
%
% INPUTS:
%   r_sc_in_flyby_frame_nb - N-body spacecraft position history in Enceladus-centric flyby frame [N x 3]
%   r_sc_in_flyby_frame_LC - Linked Conic spacecraft position history in Enceladus-centric flyby frame [N x 3]
%   r_saturn_in_flyby_frame_nb - Saturn's position history in the N-body Enceladus-centric flyby frame [N x 3]
%   full_time_out_nb       - Time vector corresponding to the states [N x 1]
%   pars                   - Parameters struct containing moon data (especially Enceladus radius)

% --- Animation State Variables ---
isPlaying = false;
playbackSpeed = 1;
playbackDirection = 1; % 1 for forward, -1 for backward
time_in_minutes = full_time_out_nb / 60;
time_in_hours = full_time_out_nb / 3600;
num_steps = length(full_time_out_nb);

% --- Find the initial index for t = -200 minutes ---
[~, initialIndex] = min(abs(time_in_minutes - (-200)));
currentIndex = initialIndex;

% --- Create Figure and Axes ---
fig_anim = figure('Name', 'Animated Enceladus-Centric Flyby', 'Color', 'w', 'Position', [150 150 1000 850]);
ax_anim = axes('Parent', fig_anim, 'Position', [0.1, 0.15, 0.8, 0.75]);
hold(ax_anim, 'on');

% --- Plot Static Elements ---
% Plot Enceladus at the center (as a transparent sphere)
radius_enceladus_display = pars.Moon.EquRad*pars.EncPlotSize; % Exaggerated radius for visualization, adjust as needed
[x_sphere, y_sphere, z_sphere] = sphere(50);
surf(ax_anim, x_sphere*radius_enceladus_display, y_sphere*radius_enceladus_display, z_sphere*radius_enceladus_display, ...
     'FaceColor', [0.678, 0.847, 0.902], 'EdgeColor', 'none', 'DisplayName', 'Enceladus', 'FaceAlpha', 0.7);

% Plot the full trajectories of the spacecraft (already in Enceladus-centric flyby frame)
plot3(ax_anim, r_sc_in_flyby_frame_nb(:,1), r_sc_in_flyby_frame_nb(:,2), r_sc_in_flyby_frame_nb(:,3), ...
      'Color', [0.8500, 0.3250, 0.0980, 0.4], 'LineWidth', 1.5, 'DisplayName', 'N-Body S/C Trajectory'); % N-body trajectory with 40% opacity
plot3(ax_anim, r_sc_in_flyby_frame_LC(:,1), r_sc_in_flyby_frame_LC(:,2), r_sc_in_flyby_frame_LC(:,3), ...
      'Color', [0, 0.4470, 0.7410], 'LineStyle', ':', 'LineWidth', 1.5, 'DisplayName', 'Linked Conic S/C Trajectory');

% --- Initialize Dynamic Plot Handles ---
h_sc_nb = plot3(ax_anim, NaN, NaN, NaN, 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.8500, 0.3250, 0.0980], 'MarkerEdgeColor', 'k', 'DisplayName', 'N-Body S/C ');
h_sc_lc = plot3(ax_anim, NaN, NaN, NaN, 's', 'MarkerSize', 6, 'MarkerFaceColor', [0, 0.4470, 0.7410], 'MarkerEdgeColor', 'k', 'DisplayName', 'Linked Conic S/C');
% NEW: Handle for Saturn's moving position
h_saturn_pos = plot3(ax_anim, NaN, NaN, NaN, 'p', 'MarkerSize', 15, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'DisplayName', 'Saturn');


% --- General Plot Formatting ---
xlabel(ax_anim, 'Flyby Frame X [km]');
ylabel(ax_anim, 'Flyby Frame Y [km]');
zlabel(ax_anim, 'Flyby Frame Z [km]');
axis(ax_anim, 'equal');
grid(ax_anim, 'on');
view(ax_anim, 135, 25);
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

% --- UI Callback Functions (Nested Functions) ---
    function onPlayPause(~,~)
        isPlaying = ~isPlaying;
        if isPlaying
            set(h_button_play, 'String', char(9616));
            if currentIndex >= num_steps || currentIndex <= 1 % If at an end, restart
                currentIndex = initialIndex;
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
        if playbackDirection == -1, playbackSpeed = 1; else, playbackSpeed = min(256, playbackSpeed * 2); end
        playbackDirection = 1;
        updateStatusText();
        if ~isPlaying
            val_from_textbox = str2double(get(h_edit_min, 'String'));
            [~, k_start] = min(abs(time_in_minutes - val_from_textbox));
            currentIndex = k_start;
            onPlayPause();
        end
    end

    function onBackward(~,~)
        if playbackDirection == 1, playbackSpeed = 1; else, playbackSpeed = min(256, playbackSpeed * 2); end
        playbackDirection = -1;
        updateStatusText();
        if ~isPlaying
            val_from_textbox = str2double(get(h_edit_min, 'String'));
            [~, k_start] = min(abs(time_in_minutes - val_from_textbox));
            currentIndex = k_start;
            onPlayPause();
        end
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

% --- Helper Functions (Nested Functions) ---
    function updateFrame(k)
        % Updates all visual elements and UI text for a given index k.
        set(h_sc_nb, 'XData', r_sc_in_flyby_frame_nb(k,1), 'YData', r_sc_in_flyby_frame_nb(k,2), 'ZData', r_sc_in_flyby_frame_nb(k,3));
        set(h_sc_lc, 'XData', r_sc_in_flyby_frame_LC(k,1), 'YData', r_sc_in_flyby_frame_LC(k,2), 'ZData', r_sc_in_flyby_frame_LC(k,3));
        % NEW: Update Saturn's position
        set(h_saturn_pos, 'XData', r_saturn_in_flyby_frame_nb(k,1), 'YData', r_saturn_in_flyby_frame_nb(k,2), 'ZData', r_saturn_in_flyby_frame_nb(k,3));
        
        set(h_title, 'String', sprintf('Enceladus-Centric Flyby | Time: %.2f minutes (%.3f hours)', time_in_minutes(k), time_in_hours(k)));
        set(h_edit_min, 'String', sprintf('%.2f', time_in_minutes(k)));
        set(h_edit_hr, 'String', sprintf('%.3f', time_in_hours(k)));
        drawnow;
    end

    function updateStatusText()
        if ~isPlaying, set(h_text_status, 'String', 'Paused'); return; end
        if playbackDirection == 1
            direction_arrow = '>>';
        else
            direction_arrow = '<<';
        end
        set(h_text_status, 'String', sprintf('Speed: %dx %s', playbackSpeed, direction_arrow));
    end

end