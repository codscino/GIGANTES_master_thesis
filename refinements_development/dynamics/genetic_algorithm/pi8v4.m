
% script for the comparison of LC and nbody trajectories for a long
% distance(64 soi)
% LC propagated from the center of Enceladus
% nbody propagated from the pericentre

clc;
clear all;
close all;

soi_multiplier = 64; % how much soi distance to stop the propagation

pars.GroundTr.t_prop  = 50*2000; % minutes, 2000 is more that 1 enceladus full revolution

perturbingBodyNaifIDs = [-2, -3, -4, 606, 605, 10, 603, 604, 601, 607, 608, 602, 5];
% perturbingBodyNaifIDs = [-2, 10, 602];
% perturbingBodyNaifIDs = [-2, 602];
% perturbingBodyNaifIDs = [10, 602];
% perturbingBodyNaifIDs = [602]; % 3body problem
% perturbingBodyNaifIDs = []; % 2 body problem

% pars.EncPlotSize = 1; %true enceladus size, usually too tiny
pars.EncPlotSize = 30;

%% ========================================================================
%  1. DEFINE INPUTS & INITIALIZE
%  ========================================================================

% --- Flyby and Central Body Parameters ---
pars.INPUTS.idCentral = 6;      % Central Body: Saturn
pars.INPUTS.idMoon    = 1;      % Enceladus flyby
pars.INPUTS.V_inf     = 4.0;    % Hyperbolic excess velocity at Enceladus [km/s]
pars.INPUTS.Flyby.min_h = 25;   % Minimum flyby altitude [km]

% --- Define Incoming and Outgoing Asymptotes (Nodes) ---
nodein =  [pars.INPUTS.V_inf, 0.15, 0];            % [km/s, rad, rad]
nodeout = [pars.INPUTS.V_inf, 0.15, deg2rad(1)];   % [km/s, rad, rad]

% --- SPICE and Time Parameters ---
kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels); 

spiceParam.frame    = 'J2000';
spiceParam.abcorr   = 'NONE';
spiceParam.observer = '699'; % Saturn Center NAIF ID
pars.INPUTS.epoch0  = date2mjd2000([2035, 1, 10, 0, 0, 0]);



actualBodyNaifIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs >= 0);
specialPerturbationIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs < 0);

% --- Load Gravitational Parameters ---
mu_central_body = getAstroConstants('Saturn', 'Mu');
[~, mu_enceladus, R_enceladus, ~] = satMoonsConstants(1); % Enceladus

% --- Dynamically build the list of gravitational parameters (mu) for actual bodies ---
mu_TBs = zeros(1, length(actualBodyNaifIDs));
for i = 1:length(actualBodyNaifIDs)
    id = actualBodyNaifIDs(i);
    switch id
        case 10,  mu_TBs(i) = getAstroConstants('Sun', 'Mu');
        case 5, mu_TBs(i) = getAstroConstants('Jupiter', 'Mu');
        case 601, [~, mu, ~, ~] = satMoonsConstants(0); mu_TBs(i) = mu; % Mimas
        case 602, [~, mu, ~, ~] = satMoonsConstants(1); mu_TBs(i) = mu; % Enceladus
        case 603, [~, mu, ~, ~] = satMoonsConstants(2); mu_TBs(i) = mu; % Tethys
        case 604, [~, mu, ~, ~] = satMoonsConstants(3); mu_TBs(i) = mu; % Dione
        case 605, [~, mu, ~, ~] = satMoonsConstants(4); mu_TBs(i) = mu; % Rhea
        case 606, [~, mu, ~, ~] = satMoonsConstants(5); mu_TBs(i) = mu; % Titan
        case 607, mu_TBs(i) = 0.374;   % Hyperion
        case 608, mu_TBs(i) = 120.4;   % Iapetus
        otherwise
            error('lc_vs_nbody_sma:UnknownNaifID', ...
                  'The NAIF ID %d is not recognized for mu calculation.', id);
    end
end

% Retrieve Saturn Parameters 
[pars.Planet.mu, pars.Planet.EquRad, pars.Planet.OrbRad, pars.Planet.hmin] = planetConstants(pars.INPUTS.idCentral);

% Retrieve Desired Moon Parameters 
if pars.INPUTS.idCentral == 6
    [pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = satMoonsConstants(pars.INPUTS.idMoon);
end
pars.Moon.Vel    = sqrt(pars.Planet.mu/pars.Moon.OrbRad);
pars.Moon.Period = 2*pi*sqrt(pars.Moon.OrbRad^3/pars.Planet.mu);
pars.Moon.HillSph = pars.Moon.OrbRad*( pars.Moon.mu/(3*(pars.Moon.mu + pars.Planet.mu)))^(1/3);

% --- Parameters required for Linked Conic Calculation ---
pars.GroundTr.npoints      = 30e3; 

% propagate to 64 SOI (stable sma)
m_enceladus = 1.08e20;
m_saturn = getAstroConstants('Saturn', 'Mass');
a_enceladus = 238020;
r_soi_enceladus = a_enceladus * (m_enceladus / m_saturn)^(2/5);

pars.INPUTS.maxPropagationDistance = r_soi_enceladus*soi_multiplier;

pars.INPUTS.Flyby.hMapping = 300;

pars.INPUTS.NAIFCentral = 699;
pars.INPUTS.NAIFMoon = 602;


%% ========================================================================
%  2. CALCULATE INITIAL PERICENTER STATE FOR PROPAGATION
%  ========================================================================

[r_enceladus_at_flyby, v_enceladus_at_flyby] = EphSS_car_spice2(602, pars.INPUTS.epoch0 , true, spiceParam);
[vvinfin, r0_sc_in, v0_sc_in, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodein(1), nodein(2), nodein(3), pars.INPUTS.epoch0, pars);
[vvinfout, r0_sc_out, v0_sc_out, vvga] = vinfAlphaCrank_to_VinfCARTClaudio(nodeout(1), nodeout(2), nodeout(3), pars.INPUTS.epoch0, pars );

rp_flyby  = pars.INPUTS.Flyby.min_h + pars.Moon.EquRad;
e_fly     = 1 + ((rp_flyby*pars.INPUTS.V_inf^2)/pars.Moon.mu);
delta_max = 2*asin(1/e_fly);
pars.delta_max = delta_max;

% b1 = -r0_sc_out./norm(r0_sc_out);
% b3 = cross(r0_sc_out, vvga)./norm(cross(r0_sc_out, vvga));
% b2 = cross(b3,b1);
% Rm = [ b1' b2' b3' ]';

Rm = buildRm(r0_sc_out,vvga); %jose method
% Rm1 = buildRm(r_enceladus_at_flyby, v_enceladus_at_flyby); %my method

[vvinfouBM, delta] = infVelBeforeDefect(vvinfin, vvinfout, delta_max, vvga);
vvinfin_bf = [Rm*vvinfin']';
vvinfouBM_bf = [Rm*vvinfouBM']';

Energy = 0.5*norm(vvinfin_bf)^2;
sma = -mu_enceladus/(2*Energy);
ecc = 1/(sin(delta/2));
rp = sma*(1 - ecc);
hhat = cross(vvinfin_bf, vvinfouBM_bf)./norm(cross(vvinfin_bf, vvinfouBM_bf));
vp = sqrt(norm(vvinfin_bf)^2 + 2*mu_enceladus/rp);

rrp_bf = rp.*(vvinfin_bf - vvinfouBM_bf)./norm(vvinfin_bf - vvinfouBM_bf);
vvp_bf = vp.*cross(hhat, rrp_bf./rp);

Rm_inv = Rm';
rrp_saturn_centric = (Rm_inv * rrp_bf')' + r_enceladus_at_flyby;
vvp_saturn_centric = (Rm_inv * vvp_bf')' + v_enceladus_at_flyby;

rrp_saturn_centric = rrp_saturn_centric(:);
vvp_saturn_centric = vvp_saturn_centric(:);

% coherence check 
[Flyby] = Flyby_BuildUp_claudio(nodein, nodeout, pars);
stin = [r0_sc_in, v0_sc_in];
stout = [r0_sc_out, v0_sc_out];
errin = stin-Flyby.State_In;
errout = stout-Flyby.State_Out;

%% ========================================================================
%  3.1 PROPAGATE THE NBODY TRAJECTORY
%  ========================================================================

duration_sec = pars.GroundTr.t_prop * 60;
time_steps = pars.GroundTr.npoints;
time_vector_fwd = linspace(eps, duration_sec, time_steps)';
time_vector_bwd = linspace(0, -duration_sec, time_steps)';

ephem_handle = @(t_sec) get_body_positions_wrapper(t_sec, pars.INPUTS.epoch0, actualBodyNaifIDs, spiceParam);

% Create a simpler handle to get ONLY Enceladus's position for the event function
ephem_enceladus = @(t) get_body_positions_wrapper(t, pars.INPUTS.epoch0, pars.INPUTS.NAIFMoon, spiceParam);

% The anonymous function for the event. It captures 'ephem_enceladus'
% and 'pars' from the current workspace.
soi_event_func = @(t, x) soiCrossingEvent(t, x, ephem_enceladus, pars.INPUTS.maxPropagationDistance);

% Propagate FORWARD until the event is triggered
[time_out_fwd, state_out_fwd, te_fwd, ~, ~] = propagateNBodyODE3(rrp_saturn_centric, vvp_saturn_centric, time_vector_fwd, mu_central_body, mu_TBs, ephem_handle, specialPerturbationIDs, soi_event_func);

% Propagate BACKWARD until the event is triggered
[time_out_bwd, state_out_bwd, te_bwd, ~, ~] = propagateNBodyODE3(rrp_saturn_centric, vvp_saturn_centric, time_vector_bwd, mu_central_body, mu_TBs, ephem_handle, specialPerturbationIDs, soi_event_func);

% --- Check and report if the events were found ---
if ~isempty(te_fwd)
    fprintf('Forward propagation stopped at t = %.2f seconds (%.2f min) when distance to Enceladus was reached.\n', te_fwd(1), te_fwd(1)/60);
else
    fprintf('Forward propagation completed without triggering the event.\n');
end
if ~isempty(te_bwd)
    fprintf('Backward propagation stopped at t = %.2f seconds (%.2f min) when distance to Enceladus was reached.\n', te_bwd(1), te_bwd(1)/60);
else
    fprintf('Backward propagation completed without triggering the event.\n');
end

time_out_bwd = flipud(time_out_bwd);
state_out_bwd = flipud(state_out_bwd);
full_time_out_nb = [time_out_bwd; time_out_fwd(2:end)]; %from 2 so the pericentre is not repeated
full_state_out_nb = [state_out_bwd; state_out_fwd(2:end, :)];
num_steps = length(full_time_out_nb);

%% ========================================================================
%  3.2 PROPAGATE THE LINKED CONICS TRAJECTORY (WITH EVENT DETECTION)
%  ========================================================================

% --- Propagate backward from the incoming node until the event ---
[timeLC_bwd, stateLC_bwd, ~, ~, ~] = propagateKeplerODE2(r0_sc_in', v0_sc_in', time_vector_bwd, mu_central_body, soi_event_func);
% [timeLC_bwd, stateLC_bwd, ~, ~, ~] = propagateKeplerODE2(rrp_saturn_centric, vvp_saturn_centric, time_vector_bwd, mu_central_body, soi_event_func);

% --- Propagate forward from the outgoing node until the event ---
[timeLC_fwd, stateLC_fwd, ~, ~, ~] = propagateKeplerODE2(r0_sc_out', v0_sc_out', time_vector_fwd, mu_central_body, soi_event_func);
% [timeLC_fwd, stateLC_fwd, ~, ~, ~] = propagateKeplerODE2(rrp_saturn_centric, vvp_saturn_centric, time_vector_fwd, mu_central_body, soi_event_func);

% Flip the backward propagation results
stateLC_bwd_flipped = flipud(stateLC_bwd);
timeLC_bwd_flipped = flipud(timeLC_bwd);

% Merge the two branches using the OUTPUTS from the ODE solver.
full_time_out_LC = [timeLC_bwd_flipped; timeLC_fwd(1:end)];
full_state_out_LC = [stateLC_bwd_flipped; stateLC_fwd(1:end, :)];



%% ========================================================================
%  3.3 RESAMPLE TRAJECTORIES TO A COMMON TIME GRID FOR COMPARISON
%  ========================================================================

% The two trajectories have different time steps. To compare them, I must
% resample one onto the time grid of the other. I will use the N-body
% time vector as the "master" reference grid.

% Use interp1 to find the LC state at the times specified by the N-body propagation.
state_LC_resampled = interp1(full_time_out_LC, full_state_out_LC, full_time_out_nb, 'pchip');

%% ========================================================================
%  4. PREPARE DATA FOR PLOTTING
%  ========================================================================

% --- Calculate the full position history of Enceladus ---
r_enc_history = zeros(num_steps, 3);
v_enc_history = zeros(num_steps, 3);
for i = 1:num_steps
    current_mjd = pars.INPUTS.epoch0 + full_time_out_nb(i) / 86400;
    [r_enc, v_enc] = EphSS_car_spice2(602, current_mjd, true, spiceParam);
    r_enc_history(i, :) = r_enc;
    v_enc_history(i, :) = v_enc;
end


%% ========================================================================
%  6. PLOT TRAJECTORIES IN ENCELADUS-CENTERED STATIC FRAME
%  ========================================================================

% Pre-allocate memory for the transformed position vectors for efficiency
r_sc_in_flyby_frame_nb = zeros(num_steps, 3);
r_sc_in_flyby_frame_LC = zeros(num_steps, 3);
r_saturn_in_flyby_frame_nb = zeros(num_steps, 3);

Rm_history = zeros(num_steps, 3, 3);


% Loop through each time step of the propagated trajectory
for i = 1:num_steps
    % --- Step 1: Translation ---
    % Get the spacecraft's position relative to Enceladus by subtracting
    % Enceladus's position from the spacecraft's position.
    % The result is still in the ICRF (J2000) frame.
    sc_pos_relative_icrf_nb = full_state_out_nb(i, 1:3) - r_enc_history(i, :);
    sc_pos_relative_icrf_LC = state_LC_resampled(i, 1:3) - r_enc_history(i, :);

    saturn_pos_relative_to_enceladus_icrf = -r_enc_history(i, :);

    % --- Step 2: Rotation ---
    % Rotate the relative position vector from the ICRF frame into the
    % flyby frame defined by the matrix Rm.
    Rm = buildRm(r_enc_history(i,:), v_enc_history(i, :));

    Rm_history(i, :, :) = Rm;

    r_sc_in_flyby_frame_nb(i, :) = (Rm * sc_pos_relative_icrf_nb')';
    r_sc_in_flyby_frame_LC(i, :) = (Rm * sc_pos_relative_icrf_LC')';

    r_saturn_in_flyby_frame_nb(i, :) = (Rm * saturn_pos_relative_to_enceladus_icrf')';
end


valid_indices_nb = vecnorm(r_sc_in_flyby_frame_nb - pars.Moon.EquRad, 2, 2)  < pars.INPUTS.maxPropagationDistance;
valid_indices_LC = vecnorm(r_sc_in_flyby_frame_LC - pars.Moon.EquRad, 2, 2)  < pars.INPUTS.maxPropagationDistance;

r_sc_in_flyby_frame_LC_valid = r_sc_in_flyby_frame_LC(valid_indices_LC,:);
r_sc_in_flyby_frame_nb_valid = r_sc_in_flyby_frame_nb(valid_indices_nb,:);

% --- Create a new figure for the 3D plot ---
figure('Name', 'Enceladus Flyby in Static Frame', 'Color', 'w');
hold on; % Hold the plot to draw multiple items

% --- Plot Enceladus as a sphere ---
% Define the radius and a light blue color for the moon
radius = pars.Moon.EquRad*pars.EncPlotSize;
light_blue_color = [0.678, 0.847, 0.902]; % light blue

% Generate the coordinates for a sphere
[x_sphere, y_sphere, z_sphere] = sphere(50); 

% Draw the sphere, scaling it by the moon's radius
surf(x_sphere*radius, y_sphere*radius, z_sphere*radius, ...
     'FaceColor', light_blue_color, ...
     'EdgeColor', 'none', ...
     'DisplayName', 'Enceladus');
alpha(0.7); % Make the sphere slightly transparent to see trajectories behind it

% --- Plot the Linked Conics Trajectory ---
plot3(r_sc_in_flyby_frame_LC(:,1), r_sc_in_flyby_frame_LC(:,2), r_sc_in_flyby_frame_LC(:,3), ...
      'Color', [0, 0.4470, 0.7410], ... % Blue color
      'LineStyle', '--', ...             % Dashed line
      'LineWidth', 2, ...
      'DisplayName', 'Linked Conic Trajectory');

% --- Plot the N-Body Trajectory ---
plot3(r_sc_in_flyby_frame_nb(:,1), r_sc_in_flyby_frame_nb(:,2), r_sc_in_flyby_frame_nb(:,3), ...
      'Color', [0.8500, 0.3250, 0.0980 0.6], ... % Orange color
      'LineWidth', 4, ...
      'DisplayName', 'N-Body Trajectory');

% --- Add plot formatting and labels ---
title('Flyby Trajectory in Enceladus-Centered Static Frame');
xlabel('Flyby Frame X [km]');
ylabel('Flyby Frame Y [km]');
zlabel('Flyby Frame Z [km]');
legend('show', 'Location', 'best'); % Display the legend
grid on;       % Add a grid
axis equal;    % Ensure the sphere looks like a sphere

% --- Set a better view and lighting ---
view(135, 25); % Set a custom viewing angle (azimuth, elevation)
camlight left;  % Add a light source from the left
lighting gouraud; % Improve the lighting on the sphere


hold off;


%% ========================================================================
%  7. DEVIATION ANALYSIS: LINKED CONICS VS. N-BODY
%  ========================================================================

% --- Calculate the Euclidean distance (norm of the difference) between the
% --- N-body and Linked Conic position vectors at each time step.
deviation_km = vecnorm(r_sc_in_flyby_frame_nb - r_sc_in_flyby_frame_LC, 2, 2);

% --- Calculate Keplerian elements for both trajectories at each time step ---
% This is needed for the semi-major axis comparison plots.
sma_nb = zeros(num_steps, 1);
sma_lc = zeros(num_steps, 1);

for i = 1:num_steps
    % Convert N-body state to Keplerian elements
    kep_nb = car2kep(full_state_out_nb(i, :), mu_central_body);
    sma_nb(i) = kep_nb(1); % Extract semi-major axis

    % Convert Linked Conic state to Keplerian elements
    kep_lc = car2kep(state_LC_resampled(i, :), mu_central_body);
    sma_lc(i) = kep_lc(1); % Extract semi-major axis
end

% Calculate the difference in semi-major axis
sma_difference_km = sma_nb - sma_lc;


% --- Create a new figure for the deviation plots ---
figure('Name', 'N-Body vs. Linked Conic Deviation Analysis', 'Color', 'w', 'Position', [100 100 1200 1200]); % Adjusted position for 3 rows

% Convert time vector from seconds to minutes for better readability
time_in_minutes = full_time_out_nb / 60;
% Calculate distance from Enceladus, normalized by its SOI radius
distance_from_enceladus_normalized = vecnorm(r_sc_in_flyby_frame_nb, 2, 2)/r_soi_enceladus;
% Mark periapsis
[min_dist, min_idx] = min(distance_from_enceladus_normalized);


% --- SUBPLOT 1: Position Deviation as a function of Time ---
ax1 = subplot(3, 2, 1); % Changed to 3 rows, 2 columns, position 1
hold(ax1, 'on');

plot(ax1, time_in_minutes, deviation_km, 'b-', 'LineWidth', 2, 'DisplayName', 'Position Deviation');
xline(ax1, 0, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Periapsis (t=0)', 'Label', 'Periapsis');

% Formatting
title(ax1, 'Position Deviation vs. Time');
xlabel(ax1, 'Time from Periapsis [minutes]');
ylabel(ax1, 'Position Deviation [km]');
legend(ax1, 'show', 'Location', 'northwest');
grid(ax1, 'on');
hold(ax1, 'off');

% --- SUBPLOT 2: Position Deviation as a function of Distance ---
ax2 = subplot(3, 2, 2); % Changed to 3 rows, 2 columns, position 2
hold(ax2, 'on');

plot(ax2, distance_from_enceladus_normalized, deviation_km, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2, 'DisplayName', 'Position Deviation');

dev_at_min_dist = deviation_km(min_idx);
plot(ax2, min_dist, dev_at_min_dist, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Periapsis');

% Add labels for Inbound and Outbound branches
text(distance_from_enceladus_normalized(1), deviation_km(1), ' Inbound', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
text(distance_from_enceladus_normalized(end), deviation_km(end), ' Outbound', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');

% Formatting
title(ax2, 'Position Deviation vs. Distance from Enceladus');
xlabel(ax2, 'Distance from Enceladus [Enceladus SOI]');
ylabel(ax2, 'Position Deviation [km]');
legend(ax2, 'show', 'Location', 'southwest');
grid(ax2, 'on');
hold(ax2, 'off');

% --- SUBPLOT 3: Semi-Major Axis Difference vs. Time ---
ax3 = subplot(3, 2, 3); % Changed to 3 rows, 2 columns, position 3
hold(ax3, 'on');

plot(ax3, time_in_minutes, sma_difference_km, 'g-', 'LineWidth', 2, 'DisplayName', 'SMA Difference');
xline(ax3, 0, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Periapsis (t=0)', 'Label', 'Periapsis');

% Formatting
title(ax3, 'Semi-Major Axis Difference vs. Time');
xlabel(ax3, 'Time from Periapsis [minutes]');
ylabel(ax3, 'SMA deviation [km]');
legend(ax3, 'show', 'Location', 'northwest');
grid(ax3, 'on');
hold(ax3, 'off');

% --- SUBPLOT 4: Semi-Major Axis Difference vs. Distance ---
ax4 = subplot(3, 2, 4); % Changed to 3 rows, 2 columns, position 4
hold(ax4, 'on');

plot(ax4, distance_from_enceladus_normalized, sma_difference_km, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2, 'DisplayName', 'SMA Difference');

sma_diff_at_min_dist = sma_difference_km(min_idx);
plot(ax4, min_dist, sma_diff_at_min_dist, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Periapsis');

% Add labels for Inbound and Outbound branches
text(distance_from_enceladus_normalized(1), sma_difference_km(1), ' Inbound', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
text(distance_from_enceladus_normalized(end), sma_difference_km(end), ' Outbound', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');

% Formatting
title(ax4, 'Semi-Major Axis Difference vs. Distance from Enceladus');
xlabel(ax4, 'Distance from Enceladus [Enceladus SOI]');
ylabel(ax4, 'SMA deviation [km]');
legend(ax4, 'show', 'Location', 'southwest');
grid(ax4, 'on');
hold(ax4, 'off');


% --- SUBPLOT 5: Raw Semi-Major Axis Values vs. Time ---
ax5 = subplot(3, 2, 5); % New subplot, 3 rows, 2 columns, position 5
hold(ax5, 'on');

plot(ax5, time_in_minutes, sma_nb, 'b-', 'LineWidth', 2, 'DisplayName', 'N-Body SMA');
plot(ax5, time_in_minutes, sma_lc, 'r--', 'LineWidth', 2, 'DisplayName', 'Linked Conics SMA');
xline(ax5, 0, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Periapsis (t=0)', 'Label', 'Periapsis');

% Formatting
title(ax5, 'Semi-Major Axis Values vs. Time');
xlabel(ax5, 'Time from Periapsis [minutes]');
ylabel(ax5, 'Semi-Major Axis [km]');
legend(ax5, 'show', 'Location', 'northwest');
grid(ax5, 'on');
hold(ax5, 'off');

% --- SUBPLOT 6: Raw Semi-Major Axis Values vs. Distance ---
ax6 = subplot(3, 2, 6); % New subplot, 3 rows, 2 columns, position 6
hold(ax6, 'on');

plot(ax6, distance_from_enceladus_normalized, sma_nb, 'b-', 'LineWidth', 2, 'DisplayName', 'N-Body SMA');
plot(ax6, distance_from_enceladus_normalized, sma_lc, 'r--', 'LineWidth', 2, 'DisplayName', 'Linked Conics SMA');

% Mark periapsis for both
sma_nb_at_min_dist = sma_nb(min_idx);
sma_lc_at_min_dist = sma_lc(min_idx);
plot(ax6, min_dist, sma_nb_at_min_dist, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'N-Body SMA at Periapsis');
plot(ax6, min_dist, sma_lc_at_min_dist, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'LC SMA at Periapsis');

% Add labels for Inbound and Outbound branches
text(distance_from_enceladus_normalized(1), sma_nb(1), ' N-Body Inbound', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold', 'Color', 'b');
text(distance_from_enceladus_normalized(1), sma_lc(1), ' LC Inbound', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'Color', 'r');

text(distance_from_enceladus_normalized(end), sma_nb(end), ' N-Body Outbound', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontWeight', 'bold', 'Color', 'b');
text(distance_from_enceladus_normalized(end), sma_lc(end), ' LC Outbound', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'Color', 'r');

% Formatting
title(ax6, 'Semi-Major Axis Values vs. Distance from Enceladus');
xlabel(ax6, 'Distance from Enceladus [Enceladus SOI]');
ylabel(ax6, 'Semi-Major Axis [km]');
legend(ax6, 'show', 'Location', 'southwest');
grid(ax6, 'on');
hold(ax6, 'off');


%% Adjust subplot layout to prevent overlap
sgtitle('N-Body vs. Linked Conic Analysis'); % Super title for the whole figure
linkaxes([ax1, ax3, ax5], 'x'); % Link x-axis of time plots
linkaxes([ax2, ax4, ax6], 'x'); % Link x-axis of distance plots


%% ========================================================================
%  8. NODEIN, NODEOUT at +- 64SOI
%  ========================================================================
nodein_nb = car2kep(full_state_out_nb(1,:), mu_central_body);
nodein_LC = car2kep(state_LC_resampled(1,:), mu_central_body);

nodeout_nb = car2kep(full_state_out_nb(end,:), mu_central_body);
nodeout_LC = car2kep(state_LC_resampled(end,:), mu_central_body);



%% ========================================================================
%  9. CREATE ANIMATED PLOTS
%  ========================================================================

% --- Animated Saturn-Centric Flyby Plot ---
createAnimatedFlybyPlot(full_state_out_nb, state_LC_resampled, r_enc_history, full_time_out_nb, pars);

% --- Animated Enceladus-Centric Flyby Plot ---
createAnimatedEnceladusCentricPlot(r_sc_in_flyby_frame_nb, r_sc_in_flyby_frame_LC, r_saturn_in_flyby_frame_nb, full_time_out_nb, pars);



%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================
function positions = get_body_positions_wrapper(t_sec, t0_mjd_start, naif_ids, sp_param)
    current_mjd = t0_mjd_start + t_sec / 86400;
    num_bodies = length(naif_ids);
    positions = zeros(3, num_bodies);
    
    for k = 1:num_bodies
        [r, ~] = EphSS_car_spice2(naif_ids(k), current_mjd, true, sp_param);
        positions(:, k) = r(:); % Ensure column vector
    end
end


function [value, isterminal, direction] = soiCrossingEvent(t, x, ephem_enceladus_handle, target_distance)
    % This function stops integration when the spacecraft's distance to
    % Enceladus equals the target distance.
    
    % 1. Get current spacecraft position from the state vector x
    r_sc = x(1:3);
    
    % 2. Get current Enceladus position using the provided handle
    r_enceladus = ephem_enceladus_handle(t);
    
    % 3. Calculate the current distance between them
    current_distance = norm(r_sc - r_enceladus);
    
    % 4. Define the event trigger. Event occurs when value = 0.
    value = current_distance - target_distance;
    
    % 5. Stop the integration when the event occurs.
    isterminal = 1;
    
    % 6. Detect the zero crossing regardless of direction (approaching or leaving).
    direction = 0;
end