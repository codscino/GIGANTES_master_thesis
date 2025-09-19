function [burn_schedule, total_deltav] = planMultiBurnStrategy(B_error_initial, max_deltav_per_burn, sensitivity_matrix)
% planMultiBurnStrategy: Plans multiple burns to achieve B-plane targeting
% when a single burn would exceed deltaV constraints
%
% INPUTS:
%   B_error_initial      - [2x1] Initial B-plane error [B_R_error; B_T_error] (km)
%   max_deltav_per_burn  - Maximum deltaV per burn (m/s)
%   sensitivity_matrix   - [2x3] Sensitivity of B-plane to velocity changes (km per km/s)
%
% OUTPUTS:
%   burn_schedule - Structure array with burn information
%   total_deltav  - Total deltaV across all burns (m/s)

fprintf('\n=== MULTI-BURN STRATEGY PLANNING ===\n');

% Convert max deltaV to km/s for calculations
max_dv_kms = max_deltav_per_burn / 1000;

% Initial error magnitude
error_magnitude = norm(B_error_initial);
fprintf('Initial B-plane error: %.3f km\n', error_magnitude);
fprintf('Max deltaV per burn: %.3f m/s\n', max_deltav_per_burn);

% Estimate total deltaV needed (unconstrained)
dv_unconstrained = pinv(sensitivity_matrix) * (-B_error_initial);
dv_unconstrained_ms = norm(dv_unconstrained) * 1000;

fprintf('Unconstrained solution would require: %.3f m/s\n', dv_unconstrained_ms);

% Determine number of burns needed
num_burns = ceil(dv_unconstrained_ms / max_deltav_per_burn);
fprintf('Number of burns needed: %d\n\n', num_burns);

% Initialize burn schedule
burn_schedule = [];
remaining_error = B_error_initial;
total_deltav = 0;

% Plan each burn
for i = 1:num_burns
    fprintf('Burn %d:\n', i);
    
    % Calculate optimal correction for remaining error
    dv_optimal = pinv(sensitivity_matrix) * (-remaining_error);
    
    % Check if we need to limit this burn
    if norm(dv_optimal) > max_dv_kms
        % Scale down to maximum allowed
        scaling = max_dv_kms / norm(dv_optimal);
        dv_applied = dv_optimal * scaling;
        fprintf('  DeltaV: %.3f m/s (constrained)\n', max_deltav_per_burn);
    else
        dv_applied = dv_optimal;
        fprintf('  DeltaV: %.3f m/s\n', norm(dv_applied) * 1000);
    end
    
    % Calculate expected B-plane change
    B_change = sensitivity_matrix * dv_applied;
    remaining_error = remaining_error + B_change;
    
    % Store burn information
    burn_schedule(i).burn_number = i;
    burn_schedule(i).deltav_vector_kms = dv_applied;
    burn_schedule(i).deltav_magnitude_ms = norm(dv_applied) * 1000;
    burn_schedule(i).B_error_after = remaining_error;
    burn_schedule(i).error_magnitude_after = norm(remaining_error);
    
    % Suggest timing (distribute burns evenly)
    burn_schedule(i).suggested_soi_distance = 64 - (i-1) * (64/num_burns);
    
    total_deltav = total_deltav + burn_schedule(i).deltav_magnitude_ms;
    
    fprintf('  B-plane error after burn: %.3f km\n', burn_schedule(i).error_magnitude_after);
    fprintf('  Suggested execution: %.1f SOI from Enceladus\n\n', burn_schedule(i).suggested_soi_distance);
    
    % Check if we've achieved sufficient accuracy
    if norm(remaining_error) < 1e-3  % 1 meter tolerance
        fprintf('Target accuracy achieved!\n');
        break;
    end
end

% Summary
fprintf('=== BURN SCHEDULE SUMMARY ===\n');
fprintf('Total burns planned: %d\n', length(burn_schedule));
fprintf('Total deltaV: %.3f m/s\n', total_deltav);
fprintf('Final B-plane error: %.3f km\n', burn_schedule(end).error_magnitude_after);

if burn_schedule(end).error_magnitude_after > 0.01  % 10 km
    fprintf('\nWARNING: Residual B-plane error is significant.\n');
    fprintf('Consider:\n');
    fprintf('  - Increasing deltaV budget\n');
    fprintf('  - Adjusting initial trajectory\n');
    fprintf('  - Accepting larger targeting errors\n');
end

% Create a burn timeline
fprintf('\n=== EXECUTION TIMELINE ===\n');
fprintf('SOI Distance | Burn # | DeltaV (m/s) | Error After (km)\n');
fprintf('-------------|--------|--------------|----------------\n');
for i = 1:length(burn_schedule)
    fprintf('%11.1f  | %6d | %12.3f | %15.3f\n', ...
        burn_schedule(i).suggested_soi_distance, ...
        burn_schedule(i).burn_number, ...
        burn_schedule(i).deltav_magnitude_ms, ...
        burn_schedule(i).error_magnitude_after);
end

end

%% Example usage function
function demonstrateMultiBurnStrategy()
    % Example: Plan multi-burn strategy for a challenging flyby
    
    % Typical B-plane error from your results
    B_error_initial = [0.681; 556.433]; % km [B_R_error; B_T_error]
    
    % Maximum deltaV per burn
    max_deltav_per_burn = 5.0; % m/s
    
    % Typical sensitivity matrix (you would get this from STM)
    % These are approximate values - replace with actual from STM
    sensitivity_matrix = [
        -0.15, 0.02, 0.01;   % dB_R/dv (km per km/s)
        0.05, -50.0, 2.0     % dB_T/dv (km per km/s)
    ];
    
    % Plan the burns
    [burn_schedule, total_deltav] = planMultiBurnStrategy(B_error_initial, max_deltav_per_burn, sensitivity_matrix);
    
    % Visualize the strategy
    figure('Name', 'Multi-Burn B-Plane Targeting Strategy');
    
    % Plot 1: B-plane error evolution
    subplot(2,1,1);
    burn_numbers = [0, [burn_schedule.burn_number]];
    errors = [norm(B_error_initial), [burn_schedule.error_magnitude_after]];
    plot(burn_numbers, errors, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Burn Number');
    ylabel('B-Plane Error (km)');
    title('B-Plane Error Reduction with Multiple Burns');
    grid on;
    
    % Add constraint line
    hold on;
    plot([0, length(burn_schedule)], [1, 1], 'r--', 'LineWidth', 1);
    legend('Error Evolution', 'Target (1 km)', 'Location', 'best');
    
    % Plot 2: DeltaV usage
    subplot(2,1,2);
    deltavs = [burn_schedule.deltav_magnitude_ms];
    bar(1:length(deltavs), deltavs);
    xlabel('Burn Number');
    ylabel('DeltaV (m/s)');
    title(sprintf('DeltaV Distribution (Total: %.2f m/s)', total_deltav));
    
    % Add constraint line
    hold on;
    plot([0, length(deltavs)+1], [max_deltav_per_burn, max_deltav_per_burn], 'r--', 'LineWidth', 2);
    legend('Burn DeltaV', 'Constraint', 'Location', 'best');
    ylim([0, max_deltav_per_burn * 1.2]);
    grid on;
end