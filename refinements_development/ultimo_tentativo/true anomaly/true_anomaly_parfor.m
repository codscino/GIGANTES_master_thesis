tic
%% STM with True Anomaly Sweep - Parallel Processing (Reduced Output)
clc;
clear all;
% close all;

%% ========================================================================
% 1. SETUP PARAMETERS
% ========================================================================

% Your existing parameter setup
pars.GroundTr.npoints = 1000;
% pars.INPUTS.perturbingBodyNaifIDs = [-2, 10];
pars.INPUTS.perturbingBodyNaifIDs = [-2,10,602, 606]; % J2, Enceladus, Sun, Titan
% pars.INPUTS.perturbingBodyNaifIDs = [];
pars.INPUTS.idCentral = 6;
pars.INPUTS.idMoon = 1;
pars.INPUTS.Flyby.min_h = 10;

% Load kernels and constants
kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels);

[pars.Planet.mu, ~, ~, ~] = planetConstants(pars.INPUTS.idCentral);
[pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = satMoonsConstants(pars.INPUTS.idMoon);
pars.Moon.HillSph = pars.Moon.OrbRad * (pars.Moon.mu / (3 * pars.Planet.mu))^(1/3);

% Other setup parameters
pars.INPUTS.NAIFCentral = 699;
pars.INPUTS.NAIFMoon = 602;
pars.INPUTS.epoch0 = date2mjd2000([2030 1 1 0 0 0]);
pars.INPUTS.V_inf = 4;

% Flyby parameters
% nodein = [4, deg2rad(45.6918), deg2rad(0)];
% nodeout = [4, deg2rad(45.6918), deg2rad(1)];

% nodein  = [4, deg2rad(8.5886), deg2rad(0)];
% nodeout = [4, deg2rad(8.6918), deg2rad(0)];

% Flyby nodes - Partial-COT 1 (O/I) - 1st Flyby
nodein = [4, deg2rad(8.6918), deg2rad(-86.9406)];
nodeout = [4, deg2rad(8.6918), deg2rad(-88.1610)];


%% ========================================================================
% 2. DEFINE SWEEP PARAMETERS
% ========================================================================

% True anomaly sweep parameters
anomaly_start = 15;
anomaly_end = 255;
anomaly_step = 25; 
anomaly_array = anomaly_start:anomaly_step:anomaly_end;
num_anomalies = length(anomaly_array);

fprintf('Starting STM sweep: %d to %d degrees (step: %d deg) - %d cases\n', ...
    anomaly_start, anomaly_end, anomaly_step, num_anomalies);

%% ========================================================================
% 3. INITIALIZE PARALLEL PROCESSING
% ========================================================================

% Start parallel pool and initialize workers
if isempty(gcp('nocreate'))
    parpool; % Start a parallel pool if one is not running
end

% Load SPICE kernels on each worker
spmd
    loadSpiceKernels(kernels);
end

%% ========================================================================
% 4. PRE-ALLOCATE STORAGE ARRAYS
% ========================================================================

% Pre-allocate storage for results
all_deltaV_magnitudes = zeros(num_anomalies, 1);
all_time_differences = zeros(num_anomalies, 1);
all_iterations = zeros(num_anomalies, 1);
all_converged_flags = zeros(num_anomalies, 1);
all_B_R_errors = zeros(num_anomalies, 1);
all_B_T_errors = zeros(num_anomalies, 1);
all_total_errors = zeros(num_anomalies, 1);
all_backward_durations = zeros(num_anomalies, 1);

%% ========================================================================
% 5. PARALLEL SWEEP LOOP
% ========================================================================

fprintf('Progress: ');
parfor anomaly_idx = 1:num_anomalies
    % Each worker gets its own copy of 'pars' and modifies it locally
    local_pars = pars;
    
    % Set the current true anomaly
    current_anomaly = anomaly_array(anomaly_idx);
    local_pars.backward_true_anomaly_deg = current_anomaly;
    
    % Show progress (note: this won't show in order due to parallel execution)
    fprintf('%.0f° ', current_anomaly);
    
    try
        % Capture console output to extract iteration count and other metrics
        original_state = warning('off', 'all'); % Suppress warnings during parallel execution
        
        % Call the STM optimization function with output suppression and capture all outputs
        [NB_in_corrected, B_achieved_final, time_diff, deltaV_mag, iterations, converged, B_R_error, B_T_error, backward_duration] = ...
            STM_tue_anomaly(local_pars, nodein, nodeout, true);
        
        % Store results directly from STM function
        all_deltaV_magnitudes(anomaly_idx) = deltaV_mag;
        all_time_differences(anomaly_idx) = time_diff;
        all_iterations(anomaly_idx) = iterations;
        all_converged_flags(anomaly_idx) = converged;
        all_B_R_errors(anomaly_idx) = B_R_error;
        all_B_T_errors(anomaly_idx) = B_T_error;
        all_total_errors(anomaly_idx) = sqrt(B_R_error^2 + B_T_error^2);
        all_backward_durations(anomaly_idx) = backward_duration;
        
        warning(original_state); % Restore warnings
        
    catch ME
        % Handle errors gracefully (silently)
        all_deltaV_magnitudes(anomaly_idx) = NaN;
        all_time_differences(anomaly_idx) = NaN;
        all_iterations(anomaly_idx) = NaN;
        all_converged_flags(anomaly_idx) = 0;
        all_B_R_errors(anomaly_idx) = NaN;
        all_B_T_errors(anomaly_idx) = NaN;
        all_total_errors(anomaly_idx) = NaN;
        all_backward_durations(anomaly_idx) = NaN;
    end
end

fprintf('\nParallel sweep completed!\n');

%% ========================================================================
% 6. PROCESS AND DISPLAY RESULTS
% ========================================================================

% Filter out failed cases
valid_indices = ~isnan(all_deltaV_magnitudes);
valid_anomalies = anomaly_array(valid_indices);
valid_deltaV = all_deltaV_magnitudes(valid_indices);
valid_time_diff = all_time_differences(valid_indices);
valid_iterations = all_iterations(valid_indices);
valid_converged = all_converged_flags(valid_indices);
valid_total_errors = all_total_errors(valid_indices);
valid_backward_durations = all_backward_durations(valid_indices);

fprintf('\nSweep Results Summary:\n');
fprintf('  Valid cases: %d / %d\n', sum(valid_indices), num_anomalies);
fprintf('  Converged cases: %d / %d\n', sum(valid_converged), sum(valid_indices));
fprintf('  DeltaV range: %.3f - %.3f m/s\n', min(valid_deltaV), max(valid_deltaV));
fprintf('  Time diff range: %.2f - %.2f minutes\n', min(valid_time_diff), max(valid_time_diff));

%% ========================================================================
% 7. CREATE VISUALIZATION
% ========================================================================

figure('Position', [100, 100, 1400, 800]);

% Subplot 1: Delta-V vs True Anomaly
subplot(2, 2, 1);
scatter(valid_anomalies, valid_deltaV, 30, valid_converged, 'filled');
colormap(gca, [1 0.2 0.2; 0.2 0.8 0.2]); % Red for not converged, Green for converged
colorbar('Ticks', [0, 1], 'TickLabels', {'Not Converged', 'Converged'});
xlabel('Backward True Anomaly [degrees]');
ylabel('Delta-V Magnitude [m/s]');
title('Delta-V vs True Anomaly');
grid on;
xlim([0, 360]);

% Subplot 2: Time Difference vs True Anomaly
subplot(2, 2, 2);
scatter(valid_anomalies, valid_time_diff, 30, valid_converged, 'filled');
colormap(gca, [1 0.2 0.2; 0.2 0.8 0.2]);
colorbar('Ticks', [0, 1], 'TickLabels', {'Not Converged', 'Converged'});
xlabel('Backward True Anomaly [degrees]');
ylabel('Time Difference [minutes]');
title('Flyby Timing vs True Anomaly');
grid on;
xlim([0, 360]);

% Subplot 3: Iterations vs True Anomaly
subplot(2, 2, 3);
scatter(valid_anomalies, valid_iterations, 30, valid_converged, 'filled');
colormap(gca, [1 0.2 0.2; 0.2 0.8 0.2]);
colorbar('Ticks', [0, 1], 'TickLabels', {'Not Converged', 'Converged'});
xlabel('Backward True Anomaly [degrees]');
ylabel('Number of Iterations');
title('Convergence Iterations vs True Anomaly');
grid on;
xlim([0, 360]);

% Subplot 4: B-plane Error vs True Anomaly
subplot(2, 2, 4);
scatter(valid_anomalies, valid_total_errors, 30, valid_converged, 'filled');
colormap(gca, [1 0.2 0.2; 0.2 0.8 0.2]);
colorbar('Ticks', [0, 1], 'TickLabels', {'Not Converged', 'Converged'});
xlabel('Backward True Anomaly [degrees]');
ylabel('Total B-plane Error [km]');
title('B-plane Error vs True Anomaly');
grid on;
xlim([0, 360]);

sgtitle(sprintf('STM Optimization Results vs Backward True Anomaly\nFlyby: V_{inf} = %.1f km/s, Min Alt = %.0f km', ...
    pars.INPUTS.V_inf, pars.INPUTS.Flyby.min_h), 'FontSize', 14, 'FontWeight', 'bold');

%% ========================================================================
% 8. FIND AND DISPLAY OPTIMAL CASES
% ========================================================================

% Find cases with minimum deltaV among converged solutions
converged_indices = valid_converged == 1;
if any(converged_indices)
    converged_deltaV = valid_deltaV(converged_indices);
    converged_anomalies = valid_anomalies(converged_indices);
    converged_time_diff = valid_time_diff(converged_indices);
    converged_iterations = valid_iterations(converged_indices);
    
    [min_deltaV, min_idx] = min(converged_deltaV);
    optimal_anomaly = converged_anomalies(min_idx);
    optimal_time_diff = converged_time_diff(min_idx);
    optimal_iterations = converged_iterations(min_idx);
    
    fprintf('\nOptimal Solution (Minimum Delta-V among converged cases):\n');
    fprintf('  True Anomaly: %.0f degrees\n', optimal_anomaly);
    fprintf('  Delta-V: %.3f m/s\n', min_deltaV);
    fprintf('  Time Difference: %.2f minutes\n', optimal_time_diff);
    fprintf('  Iterations: %d\n', optimal_iterations);
    
    % Add markers for optimal solution on plots
    for i = 1:4
        subplot(2, 2, i);
        hold on;
        switch i
            case 1
                plot(optimal_anomaly, min_deltaV, 'ko', 'MarkerSize', 12, 'LineWidth', 3);
            case 2
                plot(optimal_anomaly, optimal_time_diff, 'ko', 'MarkerSize', 12, 'LineWidth', 3);
            case 3
                plot(optimal_anomaly, optimal_iterations, 'ko', 'MarkerSize', 12, 'LineWidth', 3);
            case 4
                opt_error_idx = find(valid_anomalies == optimal_anomaly);
                plot(optimal_anomaly, valid_total_errors(opt_error_idx), 'ko', 'MarkerSize', 12, 'LineWidth', 3);
        end
        hold off;
    end
else
    fprintf('\nNo converged solutions found!\n');
end

fprintf('\n========================================\n');
fprintf('TRUE ANOMALY SWEEP COMPLETE\n');
fprintf('========================================\n');

toc