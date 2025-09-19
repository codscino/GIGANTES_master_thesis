clc;
clear all;
close all;

%% 1. Define Inputs
% ===================
% Flyby parameters
pars.INPUTS.idCentral = 6;      % Saturn
pars.INPUTS.idMoon    = 1;      % Enceladus flyby
pars.INPUTS.V_inf     = 4.0;    % km/s
pars.INPUTS.k         = 0;

% SPICE and Time parameters (you may need de440.bsp for the Sun/Jupiter)
kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels)
spiceParam.frame    = 'J2000';
spiceParam.abcorr   = 'NONE';
spiceParam.observer = '699';
t0_mjd = date2mjd2000([2035, 1, 1, 0, 0, 0]);

% Define the list of perturbations we want to study individually
% We use a cell array to hold the string 'J2' and the numeric NAIF IDs
perturbers_to_study = {'J2', 606, 605, 10, 603, 604, 601}; % J2, Titan, Rhea, Sun, Tethys, Dione, Mimas
perturber_names = {'J2 Oblateness', 'Titan', 'Rhea', 'Sun', 'Tethys', 'Dione', 'Mimas'};

% Define the range of alpha to test
alpha_range = linspace(0, 2*pi, 20);
num_alphas = length(alpha_range);

% Pre-allocate a matrix to store all our results
% Each row is an alpha value, each column is a perturber
all_errors = zeros(num_alphas, length(perturbers_to_study));

%% 2. Run the Analysis
% ====================
fprintf('Starting perturbation analysis...\n');

% Outer loop: Iterate through each perturber we want to isolate
for p_idx = 1:length(perturbers_to_study)
    
    current_perturber = perturbers_to_study{p_idx};
    fprintf('  Analyzing effect of: %s\n', perturber_names{p_idx});
    
    % Inner loop: Iterate through the alpha crank angle
    for a_idx = 1:num_alphas
        pars.INPUTS.alpha = alpha_range(a_idx);
        
        % --- Configure the flags for this specific run ---
        useJ2_flag = false;
        single_perturber_id = []; % Default to no third-body
        
        if ischar(current_perturber) && strcmpi(current_perturber, 'J2')
            % Case 1: We are testing the J2 effect
            useJ2_flag = true;
        else
            % Case 2: We are testing a single third body
            single_perturber_id = current_perturber;
        end
        
        % --- Call the main function with the configured inputs ---
        % Note: The Keplerian result (a_kep) is our baseline. It is the
        % trajectory with NO perturbations at all.
        [~, a_kep, a_nbody_perturbed, ~] = calculateApoapsisComparison( ...
            pars, spiceParam, t0_mjd, ...
            single_perturber_id, ...  % Pass either one ID or empty []
            'forward', ...
            useJ2_flag);              % Pass the J2 flag
        
        % The "error" here is the change in apoapsis due to this one perturber
        % We compare the perturbed N-body result to the unperturbed Keplerian baseline.
        perturbation_effect = a_nbody_perturbed - a_kep;
        
        all_errors(a_idx, p_idx) = perturbation_effect;
    end
end
fprintf('Analysis complete.\n');
cspice_kclear();

%% 3. Plot the Results
% ====================
figure('Name', 'Individual Perturbation Effects on Apoapsis');
hold on;
for p_idx = 1:length(perturbers_to_study)
    plot(rad2deg(alpha_range), all_errors(:, p_idx), 'LineWidth', 1.5, ...
         'DisplayName', perturber_names{p_idx});
end
hold off;

grid on;
xlabel('Crank Angle \alpha [deg]');
ylabel('Change in Apoapsis [km]');
title('Effect of Individual Perturbers on Post-Flyby Apoapsis');
legend('show', 'Location', 'best');
xlim([0 360]);