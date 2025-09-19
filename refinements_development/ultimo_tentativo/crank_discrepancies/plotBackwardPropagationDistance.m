function plotBackwardPropagationDistance(t_prop_hours)
% plotBackwardPropagationDistance: Propagates a trajectory backward in time
% for a specified duration and plots the distance to Enceladus, normalized
% by its Sphere of Influence (SOI).
%
% USAGE:
%   plotBackwardPropagationDistance(31.5); % Example for 31.5 hours
%

%% ========================================================================
%  1. SETUP & INITIAL CONDITIONS
%  ========================================================================
fprintf('--- Setting up simulation parameters ---\n');

% Clear environment (optional, but good for standalone scripts)
clc;
close all;

% Basic parameters
pars.GroundTr.npoints = 5000; % Number of points for the propagation
pars.INPUTS.epoch0 = date2mjd2000([2040 1 1 12 0 0]);
pars.INPUTS.V_inf = 4; % km/s

% SPICE parameters
spiceParam.frame = 'J2000';
spiceParam.abcorr = 'NONE';
spiceParam.observer = '699'; % Observer is Saturn
pars.INPUTS.NAIFCentral = 699;
pars.INPUTS.NAIFMoon = 602; % NAIF ID for Enceladus

% Load SPICE kernels (ensure these files are in your MATLAB path)
kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels);

% Get gravitational parameters
pars.Planet.mu = getAstroConstants('Saturn', 'Mu');
[~, pars.Moon.mu, ~, ~] = satMoonsConstants(1); % Enceladus is moon 1

% Calculate the Sphere of Influence (SOI) for Enceladus
m_enceladus = 1.08e20; % kg
m_saturn = getAstroConstants('Saturn', 'Mass'); % kg
a_enceladus = 238020; % km (semi-major axis)
r_soi_enceladus = a_enceladus * (m_enceladus / m_saturn)^(2/5);
fprintf('Enceladus SOI calculated as: %.2f km\n', r_soi_enceladus);

% Define the flyby conditions at periapsis (t=0)
% For simplicity, we assume a basic incoming trajectory
nodein = [pars.INPUTS.V_inf, 0, 0];
[~, r0_sc_in, v0_sc_in, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodein(1), nodein(2), nodein(3), pars.INPUTS.epoch0, pars);

%% ========================================================================
%  2. BACKWARD PROPAGATION
%  ========================================================================
fprintf('--- Starting backward propagation for %.2f hours ---\n', t_prop_hours);

% Define the time vector for backward propagation
duration_sec = t_prop_hours * 3600;
time_vector_bwd = linspace(0, -duration_sec, pars.GroundTr.npoints)';

% Propagate the trajectory backward using a simple 2-body model (Keplerian)
% The event function is set to [] so it propagates for the full duration.
[time_out, state_out, ~, ~, ~] = propagateKeplerODE2(r0_sc_in', v0_sc_in', time_vector_bwd, pars.Planet.mu, []);

fprintf('Propagation complete.\n');

%% ========================================================================
%  3. CALCULATE NORMALIZED DISTANCE
%  ========================================================================
fprintf('--- Calculating distance to Enceladus ---\n');

num_steps = length(time_out);
dist_to_enceladus = zeros(num_steps, 1);

% Loop through each point in the trajectory history
for i = 1:num_steps
    % Current time in MJD
    current_mjd = pars.INPUTS.epoch0 + time_out(i) / 86400;

    % Get spacecraft position from propagated state
    r_sc = state_out(i, 1:3);

    % Get Enceladus position from SPICE at the same time
    [r_enceladus, ~] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, current_mjd, true, spiceParam);

    % Calculate the distance between spacecraft and Enceladus
    dist_to_enceladus(i) = norm(r_sc - r_enceladus);
end

% Normalize the distance by the SOI radius
normalized_distance = dist_to_enceladus / r_soi_enceladus;

%% ========================================================================
%  4. PLOT THE RESULTS
%  ========================================================================
fprintf('--- Generating plot ---\n');

figure;
hold on;

% Plot the normalized distance vs. time in hours
plot(time_out / 3600, normalized_distance, 'LineWidth', 2);

% Add a horizontal line at y=1 to show the SOI boundary
yline(64, 'r--', '64 SOI Boundary', 'LineWidth', 1.5);


% Formatting
title(['Spacecraft Distance from Enceladus (Backward Propagation)']);
xlabel('Time from Periapsis (hours)');
ylabel('Distance / SOI_{Enceladus}');
grid on;
set(gca, 'FontSize', 12);
legend('S/C Trajectory', 'Location', 'best');
hold off;

fprintf('--- Analysis complete ---\n');

end