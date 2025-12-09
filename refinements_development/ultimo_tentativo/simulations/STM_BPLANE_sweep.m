%% STM Epoch Sweep Analysis with Parallel Processing
clc;
clear all;
close all;

%% ========================================================================
% 1. SETUP
% ========================================================================
soi_multiplier = 64;

% Parameter setup
pars.t_prop_hours = 25;
pars.GroundTr.npoints = 300;
pars.INPUTS.perturbingBodyNaifIDs = [-2, 602, 10]; % J2, Enceladus, Sun
pars.INPUTS.idCentral = 6;
pars.INPUTS.idMoon = 1;
pars.INPUTS.Flyby.min_h = 10;

% Load kernels for main thread
kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels);

[pars.Planet.mu, ~, ~, ~] = planetConstants(pars.INPUTS.idCentral);
[pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = satMoonsConstants(pars.INPUTS.idMoon);
pars.Moon.HillSph = pars.Moon.OrbRad * (pars.Moon.mu / (3 * pars.Planet.mu))^(1/3);

% Other setup parameters
pars.INPUTS.NAIFCentral = 699;
pars.INPUTS.NAIFMoon = 602;

% Propagate to 64 SOI
m_enceladus = 1.08e20;
m_saturn = getAstroConstants('Saturn', 'Mass');
a_enceladus = 238020;
r_soi_enceladus = a_enceladus * (m_enceladus / m_saturn)^(2/5);
pars.INPUTS.maxPropagationDistance = r_soi_enceladus*soi_multiplier;

pars.INPUTS.V_inf = 4;

% Flyby parameters
nodein = [pars.INPUTS.V_inf, 1.5, 0]; % [km/s, rad, rad]
nodeout = [pars.INPUTS.V_inf, 1.5, deg2rad(1)]; % [km/s, rad, rad]

%% ========================================================================
% 2. PARALLEL POOL SETUP WITH SPICE KERNEL LOADING
% ========================================================================

% Start parallel pool if not already running
if isempty(gcp('nocreate'))
    poolobj = parpool('local'); % Create pool with default number of workers
else
    poolobj = gcp('nocreate');
end

fprintf('Parallel pool started with %d workers\n', poolobj.NumWorkers);

% Load SPICE kernels on each worker
fprintf('Loading SPICE kernels on all workers...\n');
spmd
    loadSpiceKernels(kernels);
end
fprintf('SPICE kernels loaded successfully on all workers\n\n');

%% ========================================================================
% 3. EPOCH SWEEP SETUP
% ========================================================================
epoch_start = date2mjd2000([2030 1 1 0 0 0]);
epoch_end = date2mjd2000([2041 1 1 0 0 0]);
num_epochs = 100;

epoch_array = linspace(epoch_start, epoch_end, num_epochs);
deltaV_results = zeros(num_epochs, 1);
timeDiff_results = zeros(num_epochs, 1);
converged_flags = zeros(num_epochs, 1);
iterations_used = zeros(num_epochs, 1);
date_strings = cell(num_epochs, 1);

% Pre-compute date strings
for i = 1:num_epochs
    date_vec = mjd20002date(epoch_array(i));
    date_strings{i} = sprintf('%02d/%02d/%04d', date_vec(2), date_vec(3), date_vec(1));
end

% ========================================================================
% 4. RUN PARALLEL SWEEP WITH PROGRESS BAR
% ========================================================================
fprintf('\n========================================\n');
fprintf('STARTING PARALLEL EPOCH SWEEP ANALYSIS\n');
fprintf('========================================\n\n');
fprintf('Processing %d epochs in parallel...\n\n', num_epochs);

% Setup progress tracking
D = parallel.pool.DataQueue;
start_time = tic;

% Setup progress update function with persistent variables
afterEach(D, @(~) updateProgressBar(num_epochs, start_time));

% Run in parallel with progress tracking
parfor i = 1:num_epochs
    % Create local copy of parameters for this worker
    pars_local = pars;
    pars_local.INPUTS.epoch0 = epoch_array(i);
    
    % Initialize results for this iteration
    deltaV_local = NaN;
    timeDiff_local = NaN;
    converged_local = 0;
    iterations_local = 0;
    
    try
        % Run the STM optimization
        [~, ~, deltaV_ms, timeDiff_min, converged, iter_count] = ...
            performSTMOptimizationParallel(pars_local, nodein, nodeout, true);
        
        deltaV_local = deltaV_ms;
        timeDiff_local = timeDiff_min;
        converged_local = converged;
        iterations_local = iter_count;
        
    catch ME
        % Error occurred - results remain NaN
        % Errors stored silently to avoid cluttering progress display
    end
    
    % Store results
    deltaV_results(i) = deltaV_local;
    timeDiff_results(i) = timeDiff_local;
    converged_flags(i) = converged_local;
    iterations_used(i) = iterations_local;
    
    % Send progress update
    send(D, i);
end

elapsed_time = toc(start_time);
fprintf('\nParallel processing complete in %.2f seconds (%.2f sec/epoch average)\n', ...
    elapsed_time, elapsed_time/num_epochs);

%% ========================================================================
% 5. PLOT RESULTS
% ========================================================================
figure('Position', [100, 100, 1400, 800]);

% Filter out non-converged results
valid_idx = converged_flags == 1;

% Convert epoch array to datetime for better x-axis formatting
date_vectors = zeros(num_epochs, 6);
for i = 1:num_epochs
    date_vectors(i, :) = mjd20002date(epoch_array(i));
end
datetime_array = datetime(date_vectors);

% Plot 1: DeltaV
subplot(3,1,1);
plot(datetime_array(valid_idx), deltaV_results(valid_idx), 'b.-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
if any(~valid_idx)
    plot(datetime_array(~valid_idx), deltaV_results(~valid_idx), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
end
grid on;
ylabel('\DeltaV (m/s)', 'FontSize', 12);
title(sprintf('STM Optimization Results vs Epoch (%d workers)', poolobj.NumWorkers), ...
    'FontSize', 14, 'FontWeight', 'bold');
if any(~valid_idx)
    legend('Converged', 'Failed', 'Location', 'best');
else
    legend('Converged', 'Location', 'best');
end
xlim([datetime_array(1), datetime_array(end)]);

% Plot 2: Time Difference
subplot(3,1,2);
plot(datetime_array(valid_idx), timeDiff_results(valid_idx), 'g.-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
if any(~valid_idx)
    plot(datetime_array(~valid_idx), timeDiff_results(~valid_idx), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
end
grid on;
ylabel('Time Difference (min)', 'FontSize', 12);
if any(~valid_idx)
    legend('Converged', 'Failed', 'Location', 'best');
else
    legend('Converged', 'Location', 'best');
end
xlim([datetime_array(1), datetime_array(end)]);

% Plot 3: Iterations used
subplot(3,1,3);
plot(datetime_array(valid_idx), iterations_used(valid_idx), 'm.-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
ylabel('Iterations', 'FontSize', 12);
xlabel('Date', 'FontSize', 12);
title('Convergence Iterations', 'FontSize', 12);
xlim([datetime_array(1), datetime_array(end)]);

% Add statistics
fprintf('\n========================================\n');
fprintf('SWEEP ANALYSIS COMPLETE\n');
fprintf('========================================\n');
fprintf('Total epochs analyzed: %d\n', num_epochs);
fprintf('Converged: %d (%.1f%%)\n', sum(converged_flags), 100*sum(converged_flags)/num_epochs);
fprintf('Failed: %d (%.1f%%)\n', num_epochs - sum(converged_flags), 100*(1-sum(converged_flags)/num_epochs));

if any(valid_idx)
    fprintf('\nStatistics for converged cases:\n');
    fprintf('  DeltaV:\n');
    fprintf('    Min: %.3f m/s at %s\n', min(deltaV_results(valid_idx)), ...
        date_strings{find(deltaV_results == min(deltaV_results(valid_idx)), 1)});
    fprintf('    Max: %.3f m/s at %s\n', max(deltaV_results(valid_idx)), ...
        date_strings{find(deltaV_results == max(deltaV_results(valid_idx)), 1)});
    fprintf('    Mean: %.3f m/s, Std: %.3f m/s\n', ...
        mean(deltaV_results(valid_idx)), std(deltaV_results(valid_idx)));
    
    fprintf('  Time Difference:\n');
    fprintf('    Min: %.3f min, Max: %.3f min\n', ...
        min(timeDiff_results(valid_idx)), max(timeDiff_results(valid_idx)));
    fprintf('    Mean: %.3f min, Std: %.3f min\n', ...
        mean(abs(timeDiff_results(valid_idx))), std(abs(timeDiff_results(valid_idx))));
    
    fprintf('  Convergence:\n');
    fprintf('    Mean iterations: %.1f, Max: %d\n', ...
        mean(iterations_used(valid_idx)), max(iterations_used(valid_idx)));
end

% Save results
save('stm_epoch_sweep_parallel_results.mat', 'epoch_array', 'deltaV_results', ...
     'timeDiff_results', 'converged_flags', 'iterations_used', 'date_strings');

fprintf('\nResults saved to stm_epoch_sweep_parallel_results.mat\n');



%% =====================
% HELPER FUNCTION
% ======================
function updateProgressBar(total_items, start_time)
    persistent progress_count
    
    if isempty(progress_count)
        progress_count = 0;
    end
    
    progress_count = progress_count + 1;
    elapsed = toc(start_time);
    rate = progress_count / elapsed;
    remaining = (total_items - progress_count) / rate;
    
    % Create progress bar
    barWidth = 50;
    filled = round(barWidth * progress_count / total_items);
    barStr = ['[' repmat('=', 1, filled) repmat('-', 1, barWidth-filled) ']'];
    
    % Format times
    if elapsed < 60
        elapsedStr = sprintf('%.0fs', elapsed);
    elseif elapsed < 3600
        elapsedStr = sprintf('%.0fm %.0fs', floor(elapsed/60), mod(elapsed,60));
    else
        elapsedStr = sprintf('%.0fh %.0fm', floor(elapsed/3600), floor(mod(elapsed,3600)/60));
    end
    
    if remaining < 60
        etaStr = sprintf('%.0fs', remaining);
    elseif remaining < 3600
        etaStr = sprintf('%.0fm %.0fs', floor(remaining/60), mod(remaining,60));
    else
        etaStr = sprintf('%.0fh %.0fm', floor(remaining/3600), floor(mod(remaining,3600)/60));
    end
    
    % Update progress line
    fprintf('\r%s %3d/%d (%.1f%%) | Elapsed: %s | ETA: %s', ...
        barStr, progress_count, total_items, 100*progress_count/total_items, elapsedStr, etaStr);
    
    if progress_count == total_items
        fprintf('\n');
        progress_count = 0;  % Reset for next use
    end
end











