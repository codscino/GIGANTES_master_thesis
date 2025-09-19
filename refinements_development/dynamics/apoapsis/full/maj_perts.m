clc;
clear all;
close all;

%% 1. Define Inputs
% ===================
% Flyby parameters
pars.INPUTS.idCentral = 6;      % Saturn
pars.INPUTS.idMoon    = 1;      % Enceladus flyby

pars.INPUTS.NAIFCentral = 699;      % Saturn
pars.INPUTS.NAIFMoon    = 602;      % Enceladus flyby

%nodein Partial-COT 1 (O/I) - 1st Flyby
pars.INPUTS.V_inf     = 4.0;    % km/s
pars.INPUTS.alpha     = deg2rad(8.6918);   % pump angle [rad]
pars.INPUTS.k         = deg2rad(-86.9406); % Crank angle (fixed)

% SPICE and Time parameters
kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
spiceParam.frame    = 'J2000';
spiceParam.abcorr   = 'NONE';
spiceParam.observer = '699'; % Saturn Center

% Define the list of perturbations to study individually
perturbers_to_study = [-2, -3, -4, 606, 605, 10, 603, 604, 601, 607, 608, 5];
perturber_names = { ...
    'J2 Oblateness', 'J3 Oblateness', 'J4 Oblateness', ...
    'Titan', 'Rhea', 'Sun', 'Tethys', 'Dione', 'Mimas', ...
    'Hyperion', 'Iapetus', 'Jupiter', 'Enceladus' ...
};

% Range of DATES
jd1 = date2mjd2000([2020, 1, 1, 0, 0, 0]);
jd2 = date2mjd2000([2050, 1, 1, 0, 0, 0]);
num_dates = 1000; % Number of points in time to evaluate
date_range = linspace(jd1, jd2, num_dates);



% Pre-allocate a matrix to store all our results
all_errors = zeros(num_dates, length(perturbers_to_study));

%% 2. Setup Parallel Pool and Load Kernels on Each Worker
% ========================================================

poolobj = gcp('nocreate');
if isempty(poolobj)
    poolobj = parpool('local'); % Use default number of workers
end

% Load SPICE kernels on each worker
spmd
    loadSpiceKernels(kernels);
end

%% 3. Run the Parallel Analysis
% ==============================

fprintf('Starting parallel analysis...\n');
tic; % Start timing

% Use parfor for the outer loop (perturbers)
parfor p_idx = 1:length(perturbers_to_study)
    
    current_perturber_id = perturbers_to_study(p_idx);
    current_name = perturber_names{p_idx};
    
    % Create a temporary array for this perturber's results
    temp_errors = zeros(num_dates, 1);
    
    % Inner loop over dates (kept as regular for loop for each worker)
    for d_idx = 1:num_dates
        current_mjd = date_range(d_idx);
        
        % Main function with a single perturber at the current date
        [a_lc, a_nbody, ~] = lc_vs_nbody_sma2( ...
            pars, spiceParam, current_mjd, ...
            [current_perturber_id], 'backward');
        
        % The "error" is the change in the semi-major axis due to this one perturber
        perturbation_effect = a_nbody - a_lc;
        
        temp_errors(d_idx) = perturbation_effect;
    end
    
    % Store results in the main matrix (this assignment is handled automatically by parfor)
    all_errors(:, p_idx) = temp_errors;
    
    % Note: fprintf in parfor can be messy, but this will show progress
    fprintf('  Completed analysis for: %s (ID: %d)\n', current_name, current_perturber_id);
end

analysis_time = toc;
fprintf('Parallel analysis completed in %.2f seconds.\n', analysis_time);

% Clean up SPICE kernels on all workers
spmd
    cspice_kclear();
end
fprintf('SPICE kernels cleared on all workers.\n');

%% 4. Plot the Results
% ====================
figure('Name', 'Individual Perturbation Effects on Semi-Major Axis vs. Epoch', 'Position', [100, 100, 1000, 650]);
hold on;

% Loop through and plot each result with custom line styles
for p_idx = 1:length(perturbers_to_study)
    
    current_name = perturber_names{p_idx};
    
    % This styling logic is still useful
    line_style = '-';
    marker_style = 'none';
    marker_indices = [];
    
    if startsWith(current_name, 'J')
        line_style = '--';
    elseif strcmp(current_name, 'Hyperion') || strcmp(current_name, 'Iapetus')
        marker_style = '.';
        marker_indices = 1:10:num_dates; % Adjust marker density for date range
    end
    
    %  Plot against the date_range on the x-axis ---
    plot(date_range, all_errors(:, p_idx), ...
         'LineStyle', line_style, ...
         'Marker', marker_style, ...
         'MarkerIndices', marker_indices, ...
         'LineWidth', 1.5, ...
         'DisplayName', current_name);
end

hold off;
grid on;
box on;
ylabel('Change in Semi-Major Axis (\Deltaa) [km]');
title_str = sprintf('Effect of Individual Perturbers vs. Epoch (Fixed Alpha = %.1f deg)', rad2deg(pars.INPUTS.alpha));
title(title_str);
legend('show', 'Location', 'bestoutside');
set(gca, 'FontSize', 12);
xlim([min(date_range), max(date_range)]); % Set limits before changing ticks

%%% Create Custom Date Ticks for the X-Axis %%%

% 1. Find the start and end years from MJD range
start_date_vec = mjd20002date(jd1);
end_date_vec = mjd20002date(jd2);
start_year = start_date_vec(1);
end_year = end_date_vec(1);

% 2. Define the years timespan for the ticks
tick_interval_years = 5;
tick_years = start_year:tick_interval_years:end_year;

% 3. Convert these "tick years" back to MJD to tell MATLAB where to place them
tick_mjd_positions = zeros(size(tick_years));
for i = 1:length(tick_years)
    tick_mjd_positions(i) = date2mjd2000([tick_years(i), 1, 1, 0, 0, 0]);
end

% 4. Apply the new tick positions and set the labels to be the years
set(gca, 'XTick', tick_mjd_positions);
set(gca, 'XTickLabel', tick_years);

xlabel('Epoch [Year]'); % Update the axis label to reflect the new format

%% 5. Analyze and Plot Ranked Histogram
% ======================================

% Calculate the average of the *absolute* change for each perturber
avg_abs_effects = mean(abs(all_errors), 1);

% Sort the effects and names from largest to smallest
[sorted_effects, sort_indices] = sort(avg_abs_effects, 'descend');
sorted_names = perturber_names(sort_indices);

% Create a new figure for the bar chart
figure('Name', 'Ranked Perturbation Magnitudes', 'Position', [150, 150, 900, 600]);

% Create the bar chart using the sorted data
% Using categorical arrays for the x-axis labels is robust
x_categories = categorical(sorted_names);
% Reorder categories to match our sorting, preventing alphabetical sorting
x_categories = reordercats(x_categories, sorted_names); 

b = bar(x_categories, sorted_effects, 'FaceColor', [0.3010 0.7450 0.9330]);

% Add text labels on top of each bar
xtips = b.XEndPoints;
ytips = b.YEndPoints;
labels = string(num2str(sorted_effects(:), '%.1f'));
text(xtips, ytips, labels, 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom', 'FontSize', 10, 'FontWeight', 'bold');

% Formatting the plot
grid on;
title('Ranking of Perturbation Magnitudes', 'FontSize', 16);
ylabel('Average Absolute Change in SMA |Δa| [km]', 'FontSize', 12);
xlabel('Perturbing Body / Effect', 'FontSize', 12);
set(gca, 'FontSize', 11);
