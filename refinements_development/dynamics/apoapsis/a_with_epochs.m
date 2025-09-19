clc;
clear all;
close all;

%% 1. Define Inputs
% ------------------
% Flyby parameters
pars.INPUTS.idCentral = 6;      % Saturn
pars.INPUTS.idMoon    = 1;      % Enceladus flyby
pars.INPUTS.V_inf     = 4.0;    % km/s
pars.INPUTS.alpha     = 0.15;   % rad
pars.INPUTS.k         = 0;

% SPICE and Time parameters
kernels = {'sat441.bsp', 'naif0012.tls', 'pck00010.tpc', 'de440.bsp'};
% kernels = {'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels)
spiceParam.frame    = 'J2000';
spiceParam.abcorr   = 'NONE';
spiceParam.observer = '699'; % Saturn NAIF ID
t0_mjd = date2mjd2000([2025, 1, 1, 0, 0, 0]);

% N-Body Perturbers
% perturbingBodyNaifIDs = [602, 603, 605, 604, 601, 606, 10]; % Enceladus, Tethys, Rhea, Dione, Mimas, Titan, Sun
perturbingBodyNaifIDs = 606; %Titan

useJ2 = true;

pars.INPUTS.NAIFMoon = 602;
pars.INPUTS.NAIFCentral = 699;

%% 2. Call the Function
[a_lc, a_nbody, f_in] = lc_vs_nbody_sma(pars, spiceParam, t0_mjd, perturbingBodyNaifIDs,  'backward', useJ2);


%% 3. Results
fprintf('\n--- Comparison of Apoapsis Results ---\n');
fprintf('Linked Conics: %.2f km\n', a_lc);
fprintf('N-Body:        %.2f km\n\n', a_nbody);

err = a_lc - a_nbody;

fprintf('N-Body vs. Linked Conics Error: %.2f km\n', err);


%%  j2+titan perturbations, epoch changing
jd1 = date2mjd2000([2020 1 1 0 0 0]);
jd2 = date2mjd2000([2050 1 1 0 0 0]);

l = 20;
dates = linspace(jd1, jd2, l);

a_lc_vec = zeros(l,1);
a_nbody_vec = zeros(l,1);

err_vec = zeros(l,1);

for i = 1:l
    [a_lc_vec(i), a_nbody_vec(i), ~] = lc_vs_nbody_sma(pars,...
        spiceParam, dates(i), perturbingBodyNaifIDs,  'backward', useJ2);
    err_vec(i) = abs(a_nbody_vec(i)-a_lc_vec(i));
end

% Plot everything together
figure;
hold on;
plot(dates, err_vec, 'r', 'LineWidth', 1.5);
hold off;

xlabel('time [jd]');
ylabel('a [km]');
title('Semi-major axis discrepancies between linked conics and nbody');
legend('error', 'Location','best');
grid on;


%% single perturbation, epoch changing
% Define the list of perturbations we want to study individually
% We use a cell array to hold the string 'J2' and the numeric NAIF IDs
perturbers_to_study = {'J2', 606, 605, 10, 603, 604, 601, 607, 608}; % J2, Titan, Rhea, Sun, Tethys, Dione, Mimas
perturber_names = {'J2 Oblateness', 'Titan', 'Rhea', 'Sun', 'Tethys', 'Dione', 'Mimas', 'Hyperion', 'Iapetus'};

lp = length(perturbers_to_study);
all_errors = zeros(l, lp);

% Outer loop: Iterate through each perturber we want to isolate
for p_idx = 1:lp
    
    current_perturber = perturbers_to_study{p_idx};
    fprintf('  Analyzing effect of: %s\n', perturber_names{p_idx});
    
    % Inner loop: Iterate through the epoch
    for date_idx = 1:l

        jd = dates(date_idx);
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
        

        [a_lc, a_nbody, ~] = lc_vs_nbody_sma( ...
            pars, spiceParam, jd, ...
            single_perturber_id, ...  % Pass either one ID or empty []
            'backward', ...
            useJ2_flag);              % Pass the J2 flag
        
        all_errors(date_idx, p_idx) = a_lc-a_nbody;
    end
end
fprintf('Analysis complete.\n');

%% 3. Plot the Results
% ====================
figure('Name', 'Individual Perturbation Effects on Apoapsis');
hold on;
for p_idx = 1:lp
    plot(dates, all_errors(:, p_idx), 'LineWidth', 1.5, ...
         'DisplayName', perturber_names{p_idx});
end
hold off;

grid on;
xlabel('jd date');
ylabel('Change in sma [km]');
title('Effect of Individual Perturbers on Post-Flyby sma');
legend('show', 'Location', 'best');
