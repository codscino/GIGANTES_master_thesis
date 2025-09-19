function simple_plot_linked_conics(state_LC, r_enc_history, r_titan_history, pars)
% Creates a simple 3D plot showing only the Linked Conics trajectory
%
% INPUTS:
%   state_LC           - Linked Conics trajectory state history [N x 6]
%   r_enc_history      - Position history of Enceladus [N x 3]
%   r_titan_history    - Position history of Titan [N x 3]
%   pars               - Parameters struct containing moon and planet data

% --- Planet Size Control ---
EncPlotSize = 1; % Default to normal size

% --- Create Figure and Axes ---
fig = figure('Name', 'Linked Conics Trajectory', 'Color', 'w', 'Position', [100 100 1000 700]);
ax = axes('Parent', fig);
hold(ax, 'on');

% --- Plot Saturn at origin ---
plot3(ax, 0, 0, 0, 'o', 'MarkerSize', 15, 'MarkerFaceColor', 'y', ...
      'MarkerEdgeColor', 'k', 'DisplayName', 'Saturn');

% --- Plot moon orbits ---
plot3(ax, r_enc_history(:,1), r_enc_history(:,2), r_enc_history(:,3), ...
      'k--', 'LineWidth', 1, 'DisplayName', 'Enceladus Orbit');
plot3(ax, r_titan_history(:,1), r_titan_history(:,2), r_titan_history(:,3), ...
      'Color', [1, 0.5, 0], 'LineStyle', '--', 'LineWidth', 1, ...
      'DisplayName', 'Titan Orbit');

% --- Plot Linked Conics trajectory in blue ---
plot3(ax, state_LC(:,1), state_LC(:,2), state_LC(:,3), ...
      'b-', 'LineWidth', 2, 'DisplayName', 'Linked Conics Trajectory');

% --- Find Apocenter Index ---
% Calculate distances from Saturn for all trajectory points
distances_from_saturn = sqrt(sum(state_LC(:,1:3).^2, 2));
[~, apocenter_idx] = max(distances_from_saturn);

% --- Plot Moon Spheres ---
% Create sphere geometry
[x_sphere, y_sphere, z_sphere] = sphere(20);

% Titan radius (hardcoded from original script)
radius_titan_hardcoded = 2574; % km - Titan's actual radius

% Calculate current sphere sizes
current_radius_enceladus = pars.Moon.EquRad * EncPlotSize;
current_radius_titan = radius_titan_hardcoded * EncPlotSize/5;

% Get moon positions at apocenter moment
apocenter_enc_pos = r_enc_history(apocenter_idx, :);
apocenter_titan_pos = r_titan_history(apocenter_idx, :);

% Get spacecraft position at apocenter
apocenter_sc_pos = state_LC(apocenter_idx, 1:3);

% Plot Enceladus sphere
h_enceladus_sphere = surf(ax, x_sphere*current_radius_enceladus + apocenter_enc_pos(1), ...
                             y_sphere*current_radius_enceladus + apocenter_enc_pos(2), ...
                             z_sphere*current_radius_enceladus + apocenter_enc_pos(3), ...
                             'FaceColor', [0.678, 0.847, 0.902], 'EdgeColor', 'none', ...
                             'DisplayName', 'Enceladus', 'FaceAlpha', 0.5);

% Plot Titan sphere  
h_titan_sphere = surf(ax, x_sphere*current_radius_titan + apocenter_titan_pos(1), ...
                         y_sphere*current_radius_titan + apocenter_titan_pos(2), ...
                         z_sphere*current_radius_titan + apocenter_titan_pos(3), ...
                         'FaceColor', [0.8, 0.5, 0.2], 'EdgeColor', 'none', ...
                         'DisplayName', 'Titan', 'FaceAlpha', 0.7);

% Plot spacecraft position at apocenter
plot3(ax, apocenter_sc_pos(1), apocenter_sc_pos(2), apocenter_sc_pos(3), ...
      'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', ...
      'DisplayName', 'Spacecraft at Apocenter');

% --- Add coordinate axes ---
axis_length = 2 * pars.Moon.OrbRad;
% X-axis arrow (RED)
quiver3(ax, 0, 0, 0, axis_length, 0, 0, 'r', 'LineWidth', 1, ...
        'MaxHeadSize', 0.5, 'AutoScale', 'off', 'HandleVisibility', 'off');
text(ax, axis_length*1.1, 0, 0, 'X', 'FontSize', 12, 'FontWeight', 'bold', ...
     'Color', 'r', 'HorizontalAlignment', 'center');
% Y-axis arrow (GREEN)
quiver3(ax, 0, 0, 0, 0, axis_length, 0, 'g', 'LineWidth', 1, ...
        'MaxHeadSize', 0.5, 'AutoScale', 'off', 'HandleVisibility', 'off');
text(ax, 0, axis_length*1.1, 0, 'Y', 'FontSize', 12, 'FontWeight', 'bold', ...
     'Color', 'g', 'HorizontalAlignment', 'center');
% Z-axis arrow (BLUE)
quiver3(ax, 0, 0, 0, 0, 0, axis_length, 'Color', [0, 0, 0.8], 'LineWidth', 1, ...
        'MaxHeadSize', 0.5, 'AutoScale', 'off', 'HandleVisibility', 'off');
text(ax, 0, 0, axis_length*1.1, 'Z', 'FontSize', 12, 'FontWeight', 'bold', ...
     'Color', [0, 0, 0.8], 'HorizontalAlignment', 'center');

% --- Plot formatting ---
xlabel(ax, 'X [km]', 'FontSize', 12);
ylabel(ax, 'Y [km]', 'FontSize', 12);
zlabel(ax, 'Z [km]', 'FontSize', 12);
title(ax, 'Linked Conics Trajectory in Saturn System', 'FontSize', 14, 'FontWeight', 'bold');

axis(ax, 'equal');
grid(ax, 'on');
view(ax, 45, 30);

% Disable clipping to prevent arrows from being cut off
set(ax, 'Clipping', 'off');

% --- Set axis limits based on trajectory data ---
all_x_data = [state_LC(:,1); r_enc_history(:,1); r_titan_history(:,1)];
all_y_data = [state_LC(:,2); r_enc_history(:,2); r_titan_history(:,2)];
all_z_data = [state_LC(:,3); r_enc_history(:,3); r_titan_history(:,3)];

% Add 10% padding to limits
x_range = max(all_x_data) - min(all_x_data);
y_range = max(all_y_data) - min(all_y_data);
z_range = max(all_z_data) - min(all_z_data);

padding_x = 0.1 * x_range;
padding_y = 0.1 * y_range;
padding_z = 0.1 * z_range;

xlim(ax, [min(all_x_data) - padding_x, max(all_x_data) + padding_x]);
ylim(ax, [min(all_y_data) - padding_y, max(all_y_data) + padding_y]);
zlim(ax, [min(all_z_data) - padding_z, max(all_z_data) + padding_z*10]);

% --- Add legend ---
legend(ax, 'show', 'Location', 'best');

% --- Add UI Control for Large Planets ---
h_checkbox_large = uicontrol('Parent', fig, 'Style', 'checkbox', 'String', 'Large Planets', ...
                            'Position', [20, 20, 100, 30], 'FontSize', 10, 'Value', 0, ...
                            'Callback', @onToggleLargePlanets);

% --- Lighting ---
lighting(ax, 'gouraud');
camlight(ax);

% --- Nested Function for Large Planets Toggle ---
    function onToggleLargePlanets(source, ~)
        % Update EncPlotSize based on checkbox state
        if get(source, 'Value')
            EncPlotSize = 120; % Large planets
        else
            EncPlotSize = 1;  % Normal size
        end
        
        % Recalculate sphere sizes
        current_radius_enceladus = pars.Moon.EquRad * EncPlotSize;
        current_radius_titan = radius_titan_hardcoded * EncPlotSize/5;
        
        % Update Enceladus sphere
        set(h_enceladus_sphere, 'XData', x_sphere*current_radius_enceladus + apocenter_enc_pos(1), ...
                              'YData', y_sphere*current_radius_enceladus + apocenter_enc_pos(2), ...
                              'ZData', z_sphere*current_radius_enceladus + apocenter_enc_pos(3));
        
        % Update Titan sphere
        set(h_titan_sphere, 'XData', x_sphere*current_radius_titan + apocenter_titan_pos(1), ...
                          'YData', y_sphere*current_radius_titan + apocenter_titan_pos(2), ...
                          'ZData', z_sphere*current_radius_titan + apocenter_titan_pos(3));
    end

end