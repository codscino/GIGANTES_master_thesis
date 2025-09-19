function liveplot_epoch_sweep(all_trajectories_LC, all_times, r_enc_history, r_titan_history, epoch_array, pars, all_kepga)
% Creates an interactive plot showing spacecraft orbit shape that changes with epoch slider.
% Shows orbit shapes for the spacecraft and a red sphere at the pericenter of each trajectory.
% Also displays Keplerian orbital elements.
%
% INPUTS:
%   all_trajectories_LC    - Cell array of trajectories for each epoch [num_epochs x 1]
%   all_times             - Cell array of time vectors for each epoch [num_epochs x 1]
%   r_enc_history         - Position history of Enceladus [N x 3]
%   r_titan_history       - Position history of Titan [N x 3]
%   epoch_array           - Array of epoch values (MJD2000)
%   pars                  - Parameters struct
%   all_kepga             - Array of Keplerian elements for each epoch [num_epochs x 6]

% --- Initial Setup ---
num_epochs = length(epoch_array);
currentEpochIndex = 1;

% Convert epochs to dates for display
epoch_dates = cell(num_epochs, 1);
for i = 1:num_epochs
    date_vec = mjd20002date(epoch_array(i));
    epoch_dates{i} = sprintf('%04d-%02d-%02d', date_vec(1), date_vec(2), date_vec(3));
end

% --- Create Figure and Axes ---
fig = figure('Name', 'Epoch-Dependent Spacecraft Orbit Visualization', ...
            'Color', 'w', 'Position', [100 100 1400 900]);
ax = axes('Parent', fig, 'Position', [0.05, 0.25, 0.65, 0.7]);
hold(ax, 'on');

% --- Plot Static Elements ---
% Saturn at origin
plot3(ax, 0, 0, 0, 'o', 'MarkerSize', 15, 'MarkerFaceColor', [1 0.8 0], ...
      'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', 'Saturn');

% Coordinate axes
axis_length = 2 * pars.Moon.OrbRad;
% X-axis (RED)
quiver3(ax, 0, 0, 0, axis_length, 0, 0, 'r', 'LineWidth', 1.5, ...
        'MaxHeadSize', 0.3, 'AutoScale', 'off');
text(ax, axis_length*1.1, 0, 0, 'X', 'FontSize', 14, 'FontWeight', 'bold', ...
     'Color', 'r', 'HorizontalAlignment', 'center');
% Y-axis (GREEN)
quiver3(ax, 0, 0, 0, 0, axis_length, 0, 'g', 'LineWidth', 1.5, ...
        'MaxHeadSize', 0.3, 'AutoScale', 'off');
text(ax, 0, axis_length*1.1, 0, 'Y', 'FontSize', 14, 'FontWeight', 'bold', ...
     'Color', 'g', 'HorizontalAlignment', 'center');
% Z-axis (BLUE)
quiver3(ax, 0, 0, 0, 0, 0, axis_length*0.5, 'b', 'LineWidth', 1.5, ...
        'MaxHeadSize', 0.3, 'AutoScale', 'off');
text(ax, 0, 0, axis_length*0.5*1.1, 'Z', 'FontSize', 14, 'FontWeight', 'bold', ...
     'Color', 'b', 'HorizontalAlignment', 'center');

% Plot Enceladus orbit (static)
plot3(ax, r_enc_history(:,1), r_enc_history(:,2), r_enc_history(:,3), ...
      'Color', [0.4 0.7 0.9], 'LineWidth', 0.5, 'LineStyle', '-', ...
      'DisplayName', 'Enceladus Orbit');

% Plot Titan orbit (static)
plot3(ax, r_titan_history(:,1), r_titan_history(:,2), r_titan_history(:,3), ...
      'Color', [0.8 0.4 0.2], 'LineWidth', 0.5, 'LineStyle', '-', ...
      'DisplayName', 'Titan Orbit');

% --- Initialize Dynamic Plot Handles ---
% Spacecraft Trajectory
h_sc_trajectory = plot3(ax, NaN, NaN, NaN, ...
                       'Color', [0 1 0], 'LineWidth', 0.7, ...
                       'DisplayName', 'Spacecraft Trajectory');
                       
% Pericenter Sphere
h_pericenter = plot3(ax, NaN, NaN, NaN, 'o', ...
                     'MarkerSize', 8, 'MarkerFaceColor', 'r', ...
                     'MarkerEdgeColor', 'k', 'DisplayName', 'Pericenter');

% Eccentricity Vector
h_ecc_vector = plot3(ax, [0 NaN], [0 NaN], [0 NaN], ...
                     'Color', [0 0 0], 'LineWidth', 2, ...
                     'DisplayName', 'Eccentricity Vector (scaled)');

% --- Plot Formatting ---
xlabel(ax, 'X [km]', 'FontSize', 12);
ylabel(ax, 'Y [km]', 'FontSize', 12);
zlabel(ax, 'Z [km]', 'FontSize', 12);
title(ax, '', 'FontSize', 14, 'FontWeight', 'bold');
legend(ax, 'show', 'Location', 'northeast');
grid(ax, 'on');
axis(ax, 'equal');
view(ax, 45, 30);

% Set axis limits to show the region of interest
max_orbit_radius = max([max(vecnorm(r_titan_history, 2, 2)), ...
                        max(vecnorm(r_enc_history, 2, 2))]) * 1.5;
xlim(ax, [-max_orbit_radius, max_orbit_radius]);
ylim(ax, [-max_orbit_radius, max_orbit_radius]);
zlim(ax, [-max_orbit_radius/4, max_orbit_radius/4]);

% --- Create Keplerian Elements Display Panel ---
kepPanel = uipanel('Parent', fig, 'Title', 'Enceladus Keplerian Elements', ...
                   'Position', [0.72, 0.45, 0.26, 0.5], ...
                   'FontSize', 12, 'FontWeight', 'bold');

% Keplerian elements labels and value displays
kep_labels = {'a (km):', 'e:', 'i (deg):', 'ω (deg):', 'Ω (deg):', 'θ (deg):'};
kep_handles = cell(6, 1);

for i = 1:6
    % Label
    uicontrol('Parent', kepPanel, 'Style', 'text', ...
             'String', kep_labels{i}, ...
             'Position', [10, 340 - i*50, 80, 25], ...
             'FontSize', 11, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'left');
    
    % Value display
    kep_handles{i} = uicontrol('Parent', kepPanel, 'Style', 'text', ...
                              'String', '0.000', ...
                              'Position', [100, 340 - i*50, 120, 25], ...
                              'FontSize', 11, ...
                              'HorizontalAlignment', 'left', ...
                              'BackgroundColor', [0.95 0.95 0.95]);
end

% --- Create UI Controls ---
% Epoch slider
h_slider = uicontrol('Parent', fig, 'Style', 'slider', ...
                    'Position', [150, 80, 600, 30], ...
                    'Min', 1, 'Max', num_epochs, 'Value', 1, ...
                    'SliderStep', [1/(num_epochs-1), 10/(num_epochs-1)], ...
                    'Callback', @onSliderChange);

% Epoch display text
h_text_epoch = uicontrol('Parent', fig, 'Style', 'text', ...
                        'Position', [350, 130, 200, 100], ...
                        'FontSize', 12, 'FontWeight', 'bold', ...
                        'BackgroundColor', 'w');

% Previous/Next buttons
h_btn_prev = uicontrol('Parent', fig, 'Style', 'pushbutton', ...
                      'String', '< Previous', ...
                      'Position', [50, 80, 80, 30], ...
                      'FontSize', 10, ...
                      'Callback', @onPrevious);

h_btn_next = uicontrol('Parent', fig, 'Style', 'pushbutton', ...
                      'String', 'Next >', ...
                      'Position', [770, 80, 80, 30], ...
                      'FontSize', 10, ...
                      'Callback', @onNext);

% Jump to epoch controls
uicontrol('Parent', fig, 'Style', 'text', ...
         'String', 'Jump to Year:', ...
         'Position', [900, 110, 100, 20], ...
         'HorizontalAlignment', 'right', ...
         'FontSize', 10);

h_edit_year = uicontrol('Parent', fig, 'Style', 'edit', ...
                       'Position', [1010, 105, 60, 25], ...
                       'FontSize', 10, ...
                       'Callback', @onEditYear);

% Animation controls
h_btn_animate = uicontrol('Parent', fig, 'Style', 'pushbutton', ...
                         'String', 'Animate', ...
                         'Position', [450, 40, 100, 30], ...
                         'FontSize', 11, ...
                         'Callback', @onAnimate);

% --- Initialize Display ---
updateDisplay(currentEpochIndex);

% --- Callback Functions (Nested) ---
    function onSliderChange(source, ~)
        newIndex = round(get(source, 'Value'));
        if newIndex ~= currentEpochIndex
            currentEpochIndex = newIndex;
            updateDisplay(currentEpochIndex);
        end
    end

    function onPrevious(~, ~)
        if currentEpochIndex > 1
            currentEpochIndex = currentEpochIndex - 1;
            set(h_slider, 'Value', currentEpochIndex);
            updateDisplay(currentEpochIndex);
        end
    end

    function onNext(~, ~)
        if currentEpochIndex < num_epochs
            currentEpochIndex = currentEpochIndex + 1;
            set(h_slider, 'Value', currentEpochIndex);
            updateDisplay(currentEpochIndex);
        end
    end

    function onEditYear(source, ~)
        year = str2double(get(source, 'String'));
        if ~isnan(year) && year >= 2030 && year <= 2041
            % Find closest epoch to the entered year
            target_mjd = date2mjd2000([year 1 1 0 0 0]);
            [~, newIndex] = min(abs(epoch_array - target_mjd));
            currentEpochIndex = newIndex;
            set(h_slider, 'Value', currentEpochIndex);
            updateDisplay(currentEpochIndex);
        else
            % Reset to current value if invalid
            date_vec = mjd20002date(epoch_array(currentEpochIndex));
            set(source, 'String', num2str(date_vec(1)));
        end
    end

    function onAnimate(~, ~)
        % Animate through all epochs
        set(h_btn_animate, 'Enable', 'off', 'String', 'Animating...');
        
        for idx = 1:num_epochs
            if ~isvalid(fig)
                break;
            end
            currentEpochIndex = idx;
            set(h_slider, 'Value', currentEpochIndex);
            updateDisplay(currentEpochIndex);
            pause(0.05); % Adjust speed as needed
        end
        
        if isvalid(fig)
            set(h_btn_animate, 'Enable', 'on', 'String', 'Animate');
        end
    end

    function updateDisplay(idx)
        % Update the spacecraft trajectory
        trajectory = all_trajectories_LC{idx};
        set(h_sc_trajectory, 'XData', trajectory(:,1), ...
                           'YData', trajectory(:,2), ...
                           'ZData', trajectory(:,3));
        
        % Find the middle point of the trajectory to mark the pericenter
        mid_index = round(size(trajectory, 1) / 2);
        pericenter_pos = trajectory(mid_index, 1:3);
        
        % Update the pericenter sphere's position
        set(h_pericenter, 'XData', pericenter_pos(1), ...
                        'YData', pericenter_pos(2), ...
                        'ZData', pericenter_pos(3));
        
        % Update Keplerian elements display
        kepga = all_kepga(idx, :);
        
        % Extract Keplerian elements
        a = kepga(1);     % semi-major axis
        e = kepga(2);     % eccentricity
        i = kepga(3);     % inclination
        omega = kepga(4); % argument of periapsis
        Omega = kepga(5); % longitude of ascending node
        theta = kepga(6); % true anomaly
        
        % Calculate eccentricity vector in inertial frame
        % The eccentricity vector points toward periapsis (θ = 0)
        ecc_magnitude = e * a *10000; % Scale by semi-major axis for visibility
        ecc_x = ecc_magnitude * (cos(Omega)*cos(omega) - sin(Omega)*sin(omega)*cos(i));
        ecc_y = ecc_magnitude * (sin(Omega)*cos(omega) + cos(Omega)*sin(omega)*cos(i));
        ecc_z = ecc_magnitude * (sin(omega)*sin(i));
        
        % Update eccentricity vector plot
        set(h_ecc_vector, 'XData', [0 ecc_x], ...
                         'YData', [0 ecc_y], ...
                         'ZData', [0 ecc_z]);
        
        % Format and display each Keplerian element
        set(kep_handles{1}, 'String', sprintf('%.3f', kepga(1)));           % a (km)
        set(kep_handles{2}, 'String', sprintf('%.6f', kepga(2)));           % e
        set(kep_handles{3}, 'String', sprintf('%.3f', rad2deg(kepga(3))));  % i (deg)
        set(kep_handles{4}, 'String', sprintf('%.3f', rad2deg(kepga(4))));  % ω (deg)
        set(kep_handles{5}, 'String', sprintf('%.3f', rad2deg(kepga(5))));  % Ω (deg)
        set(kep_handles{6}, 'String', sprintf('%.3f', rad2deg(kepga(6))));  % θ (deg)
        
        % Update title and text displays
        epoch_date = epoch_dates{idx};
        
        title_str = sprintf('Spacecraft Orbit at Flyby Epoch: %s', epoch_date);
        title(ax, title_str);
        
        epoch_info = sprintf('Epoch %d/%d\n%s\nMJD2000: %.2f', ...
                           idx, num_epochs, epoch_date, epoch_array(idx));
        set(h_text_epoch, 'String', epoch_info);
        
        % Update year edit box
        date_vec = mjd20002date(epoch_array(idx));
        set(h_edit_year, 'String', num2str(date_vec(1)));
        
        drawnow;
    end

end