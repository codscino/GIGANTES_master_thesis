function plot_static_three_trajectories(state_lc_merged, state_stm_merged, state_orig_merged, r_enc_history, r_titan_history, time_lc_merged, backward_duration, pars)
% Creates a static 3D plot showing three trajectories and planet orbits
% with marked positions at flyby time (t=0) and earlier epoch
%
% INPUTS:
%   state_lc_merged    - Linked Conics trajectory state history [N x 6]
%   state_stm_merged   - STM-corrected N-body trajectory state history [N x 6]
%   state_orig_merged  - Original N-body trajectory state history [N x 6]
%   r_enc_history      - Position history of Enceladus [N x 3]
%   r_titan_history    - Position history of Titan [N x 3]
%   time_lc_merged     - Time vector in seconds [N x 1]
%   backward_duration  - Duration in hours for the earlier epoch
%   pars               - Parameters struct containing moon and planet data

% Convert time to hours for easier interpretation
time_hours = time_lc_merged / 3600;

% Find indices for the two key times
[~, idx_flyby] = min(abs(time_hours)); % t = 0 (flyby time)
[~, idx_early] = min(abs(time_hours + backward_duration)); % earlier epoch

fprintf('Flyby time index: %d (t = %.3f hours)\n', idx_flyby, time_hours(idx_flyby));
fprintf('Early epoch index: %d (t = %.3f hours)\n', idx_early, time_hours(idx_early));

% Create figure
fig = figure('Name', 'Static Three-Trajectory Comparison', 'Color', 'w', 'Position', [100 100 1400 1000]);
ax = axes('Parent', fig, 'Position', [0.08, 0.12, 0.85, 0.8]);
hold(ax, 'on');

% Plot Saturn at origin
plot3(ax, 0, 0, 0, 'o', 'MarkerSize', 25, 'MarkerFaceColor', 'y', ...
      'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', 'Saturn');

% Add coordinate axes
axis_length = 1.5 * pars.Moon.OrbRad;
quiver3(ax, 0, 0, 0, axis_length, 0, 0, 'r', 'LineWidth', 2, ...
        'MaxHeadSize', 0.3, 'AutoScale', 'off');
text(ax, axis_length*1.2, 0, 0, 'X', 'FontSize', 18, 'FontWeight', 'bold', ...
     'Color', 'r', 'HorizontalAlignment', 'center');
quiver3(ax, 0, 0, 0, 0, axis_length, 0, 'g', 'LineWidth', 2, ...
        'MaxHeadSize', 0.3, 'AutoScale', 'off');
text(ax, 0, axis_length*1.2, 0, 'Y', 'FontSize', 18, 'FontWeight', 'bold', ...
     'Color', 'g', 'HorizontalAlignment', 'center');
quiver3(ax, 0, 0, 0, 0, 0, axis_length, 'b', 'LineWidth', 2, ...
        'MaxHeadSize', 0.3, 'AutoScale', 'off');
text(ax, 0, 0, axis_length*1.2, 'Z', 'FontSize', 18, 'FontWeight', 'bold', ...
     'Color', 'b', 'HorizontalAlignment', 'center');

% Plot planet orbits
plot3(ax, r_enc_history(:,1), r_enc_history(:,2), r_enc_history(:,3), ...
      'k--', 'LineWidth', 1.5, 'DisplayName', 'Enceladus Orbit');
plot3(ax, r_titan_history(:,1), r_titan_history(:,2), r_titan_history(:,3), ...
      'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Titan Orbit');

% Plot the three spacecraft trajectories
plot3(ax, state_lc_merged(:,1), state_lc_merged(:,2), state_lc_merged(:,3), ...
      'Color', [0, 0.4470, 0.7410], 'LineStyle', '-', 'LineWidth', 3, 'DisplayName', 'Linked Conics');
plot3(ax, state_stm_merged(:,1), state_stm_merged(:,2), state_stm_merged(:,3), ...
      'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3, 'DisplayName', 'STM-Corrected N-Body');
plot3(ax, state_orig_merged(:,1), state_orig_merged(:,2), state_orig_merged(:,3), ...
      'Color', [0.9290, 0.6940, 0.1250], 'LineStyle', '-', 'LineWidth', 3, 'DisplayName', 'Original N-Body');

% Plot Enceladus and Titan positions at both epochs
% At flyby time (t=0)
enceladus_flyby = plot3(ax, r_enc_history(idx_flyby,1), r_enc_history(idx_flyby,2), r_enc_history(idx_flyby,3), ...
                       'o', 'MarkerSize', 12, 'MarkerFaceColor', [0.678, 0.847, 0.902], ...
                       'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', 'Enceladus @ Flyby');
titan_flyby = plot3(ax, r_titan_history(idx_flyby,1), r_titan_history(idx_flyby,2), r_titan_history(idx_flyby,3), ...
                   'o', 'MarkerSize', 15, 'MarkerFaceColor', [0.8, 0.5, 0.2], ...
                   'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', 'Titan @ Flyby');

% At early epoch
enceladus_early = plot3(ax, r_enc_history(idx_early,1), r_enc_history(idx_early,2), r_enc_history(idx_early,3), ...
                       's', 'MarkerSize', 10, 'MarkerFaceColor', [0.678, 0.847, 0.902], ...
                       'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', 'Enceladus @ Early Epoch');
titan_early = plot3(ax, r_titan_history(idx_early,1), r_titan_history(idx_early,2), r_titan_history(idx_early,3), ...
                   's', 'MarkerSize', 13, 'MarkerFaceColor', [0.8, 0.5, 0.2], ...
                   'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', 'Titan @ Early Epoch');

% Plot spacecraft positions at both epochs
% At flyby time (t=0) - circles
plot3(ax, state_lc_merged(idx_flyby,1), state_lc_merged(idx_flyby,2), state_lc_merged(idx_flyby,3), ...
      'o', 'MarkerSize', 10, 'MarkerFaceColor', [0, 0.4470, 0.7410], ...
      'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', 'LC S/C @ Flyby');
plot3(ax, state_stm_merged(idx_flyby,1), state_stm_merged(idx_flyby,2), state_stm_merged(idx_flyby,3), ...
      'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.8500, 0.3250, 0.0980], ...
      'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', 'STM S/C @ Flyby');
plot3(ax, state_orig_merged(idx_flyby,1), state_orig_merged(idx_flyby,2), state_orig_merged(idx_flyby,3), ...
      'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.9290, 0.6940, 0.1250], ...
      'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', 'Orig S/C @ Flyby');

% At early epoch - squares
plot3(ax, state_lc_merged(idx_early,1), state_lc_merged(idx_early,2), state_lc_merged(idx_early,3), ...
      's', 'MarkerSize', 8, 'MarkerFaceColor', [0, 0.4470, 0.7410], ...
      'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', 'LC S/C @ Early');
plot3(ax, state_stm_merged(idx_early,1), state_stm_merged(idx_early,2), state_stm_merged(idx_early,3), ...
      's', 'MarkerSize', 8, 'MarkerFaceColor', [0.8500, 0.3250, 0.0980], ...
      'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', 'STM S/C @ Early');
plot3(ax, state_orig_merged(idx_early,1), state_orig_merged(idx_early,2), state_orig_merged(idx_early,3), ...
      's', 'MarkerSize', 8, 'MarkerFaceColor', [0.9290, 0.6940, 0.1250], ...
      'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', 'Orig S/C @ Early');

% Calculate and display trajectory separations at both epochs
fprintf('\nTrajectory Separations:\n');
fprintf('At Flyby Time (t = %.3f hrs):\n', time_hours(idx_flyby));
dist_lc_stm_flyby = norm(state_lc_merged(idx_flyby,1:3) - state_stm_merged(idx_flyby,1:3));
dist_lc_orig_flyby = norm(state_lc_merged(idx_flyby,1:3) - state_orig_merged(idx_flyby,1:3));
dist_stm_orig_flyby = norm(state_stm_merged(idx_flyby,1:3) - state_orig_merged(idx_flyby,1:3));
fprintf('  LC vs STM: %.3f km\n', dist_lc_stm_flyby);
fprintf('  LC vs Orig: %.3f km\n', dist_lc_orig_flyby);
fprintf('  STM vs Orig: %.3f km\n', dist_stm_orig_flyby);

fprintf('At Early Epoch (t = %.3f hrs):\n', time_hours(idx_early));
dist_lc_stm_early = norm(state_lc_merged(idx_early,1:3) - state_stm_merged(idx_early,1:3));
dist_lc_orig_early = norm(state_lc_merged(idx_early,1:3) - state_orig_merged(idx_early,1:3));
dist_stm_orig_early = norm(state_stm_merged(idx_early,1:3) - state_orig_merged(idx_early,1:3));
fprintf('  LC vs STM: %.3f km\n', dist_lc_stm_early);
fprintf('  LC vs Orig: %.3f km\n', dist_lc_orig_early);
fprintf('  STM vs Orig: %.3f km\n', dist_stm_orig_early);

% Add text annotations for the key times
text(ax, state_stm_merged(idx_flyby,1)*1.1, state_stm_merged(idx_flyby,2)*1.1, state_stm_merged(idx_flyby,3)*1.1, ...
     sprintf('Flyby\nt = %.1f hrs', time_hours(idx_flyby)), 'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', 'w', 'EdgeColor', 'k', 'HorizontalAlignment', 'center');

text(ax, state_stm_merged(idx_early,1)*1.1, state_stm_merged(idx_early,2)*1.1, state_stm_merged(idx_early,3)*1.1, ...
     sprintf('Early Epoch\nt = %.1f hrs', time_hours(idx_early)), 'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', 'w', 'EdgeColor', 'k', 'HorizontalAlignment', 'center');

% Format the plot
xlabel(ax, 'X [km]', 'FontSize', 14, 'FontWeight', 'bold');
ylabel(ax, 'Y [km]', 'FontSize', 14, 'FontWeight', 'bold');
zlabel(ax, 'Z [km]', 'FontSize', 14, 'FontWeight', 'bold');
title(ax, sprintf('Three-Trajectory Comparison\nFlyby at t=%.1f hrs, Early Epoch at t=%.1f hrs', ...
                  time_hours(idx_flyby), time_hours(idx_early)), 'FontSize', 16, 'FontWeight', 'bold');

axis(ax, 'equal');
grid(ax, 'on');
view(ax, 45, 30);
lighting gouraud;
camlight;

% Create a more organized legend
legend(ax, 'show', 'Location', 'eastoutside', 'FontSize', 10);

% Set axis limits for better visualization
axis_limits = 1.2 * max([max(abs(state_stm_merged(:,1:3)), [], 'all'), ...
                        max(abs(r_enc_history), [], 'all'), ...
                        max(abs(r_titan_history), [], 'all')]);
xlim(ax, [-axis_limits, axis_limits]);
ylim(ax, [-axis_limits, axis_limits]);
zlim(ax, [-axis_limits, axis_limits]);

% Add a subtitle with trajectory information
subtitle_text = sprintf('Circles: Positions at Flyby | Squares: Positions at Early Epoch\nBlue: Linked Conics | Orange: STM-Corrected | Yellow: Original N-Body');
annotation('textbox', [0.08, 0.02, 0.85, 0.08], 'String', subtitle_text, ...
          'FontSize', 11, 'HorizontalAlignment', 'center', ...
          'BackgroundColor', 'w', 'EdgeColor', 'k');

fprintf('\nPlot created successfully!\n');
end