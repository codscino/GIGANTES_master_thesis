clc;
clear all;
close all;

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
pars.INPUTS.epoch0  = date2mjd2000([2035, 1, 1, 0, 0, 0]);

% perturbingBodyNaifIDs = [-2, -3, -4, 606, 605, 10, 603, 604, 601, 607, 608, 602, 5];
perturbingBodyNaifIDs = [602];

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
pars.GroundTr.t_prop       = 2000;    % min, more that one enceladus full revolution


% propagate to 64 SOI (stable sma)
m_enceladus = 1.08e20;
m_saturn = getAstroConstants('Saturn', 'Mass');
a_enceladus = 238020;
r_soi_enceladus = a_enceladus * (m_enceladus / m_saturn)^(2/5);
pars.INPUTS.Flyby.hMapping = r_soi_enceladus*3064;


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
% b3 = cross(r0_sc_out, vvga)./norm(cross(r_sc_out, vvga));
% b2 = cross(b3,b1);
% Rm = [ b1' b2' b3' ]';

Rm = buildRm(r0_sc_out,vvga);

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
soi_event_func = @(t, x) soiCrossingEvent(t, x, ephem_enceladus, pars.INPUTS.Flyby.hMapping);

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

% --- Propagate forward from the outgoing node until the event ---
[timeLC_fwd, stateLC_fwd, ~, ~, ~] = propagateKeplerODE2(r0_sc_out', v0_sc_out', time_vector_fwd, mu_central_body, soi_event_func);

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
state_LC_resampled = interp1(full_time_out_LC, full_state_out_LC, full_time_out_nb, 'spline');

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

% Pre-allocate memory for the transformed state vectors for efficiency
r_sc_in_flyby_frame_nb = zeros(num_steps, 3);
v_sc_in_flyby_frame_nb = zeros(num_steps, 3); % Added for velocity
r_sc_in_flyby_frame_LC = zeros(num_steps, 3);
r_saturn_in_flyby_frame_nb = zeros(num_steps, 3);

Rm_nb_history = zeros(num_steps, 3, 3);
Rm_LC_history = zeros(num_steps, 3, 3);

% Loop through each time step of the propagated trajectory
for i = 1:num_steps
    % --- Step 1: Translation (Position & Velocity) ---
    % Get the s/c state relative to Enceladus in the ICRF (J2000) frame.
    sc_pos_relative_icrf_nb = full_state_out_nb(i, 1:3) - r_enc_history(i, :);
    sc_vel_relative_icrf_nb = full_state_out_nb(i, 4:6) - v_enc_history(i, :); % Added
    sc_pos_relative_icrf_LC = state_LC_resampled(i, 1:3) - r_enc_history(i, :);

    saturn_pos_relative_to_enceladus_icrf = -r_enc_history(i, :);

    % --- Step 2: Rotation ---
    % Build rotation matrix and rotate state vectors into the flyby frame.
    Rm_nb = buildRm(full_state_out_nb(i, 1:3), v_enc_history(i, :));
    Rm_LC = buildRm(state_LC_resampled(i, 1:3), v_enc_history(i, :));

    % Store the calculated rotation matrices
    Rm_nb_history(i, :, :) = Rm_nb;
    Rm_LC_history(i, :, :) = Rm_LC;

    r_sc_in_flyby_frame_nb(i, :) = (Rm_nb * sc_pos_relative_icrf_nb')';
    v_sc_in_flyby_frame_nb(i, :) = (Rm_nb * sc_vel_relative_icrf_nb')'; % Added
    r_sc_in_flyby_frame_LC(i, :) = (Rm_LC * sc_pos_relative_icrf_LC')';

    r_saturn_in_flyby_frame_nb(i, :) = (Rm_nb * saturn_pos_relative_to_enceladus_icrf')';
end


valid_indices_nb = vecnorm(r_sc_in_flyby_frame_nb - pars.Moon.EquRad, 2, 2)  < pars.INPUTS.Flyby.hMapping;
valid_indices_LC = vecnorm(r_sc_in_flyby_frame_LC - pars.Moon.EquRad, 2, 2)  < pars.INPUTS.Flyby.hMapping;

r_sc_in_flyby_frame_LC_valid = r_sc_in_flyby_frame_LC(valid_indices_LC,:);
r_sc_in_flyby_frame_nb_valid = r_sc_in_flyby_frame_nb(valid_indices_nb,:);

% --- Create a new figure for the 3D plot ---
figure('Name', 'Enceladus Flyby in Static Frame', 'Color', 'w');
hold on; % Hold the plot to draw multiple items

% --- Plot Enceladus as a sphere ---
% Define the radius and a light blue color for the moon
radius = pars.Moon.EquRad*30;
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
figure('Name', 'N-Body vs. Linked Conic Deviation Analysis', 'Color', 'w', 'Position', [100 100 1200 900]);

% --- SUBPLOT 1: Position Deviation as a function of Time ---
ax1 = subplot(2, 2, 1);
hold(ax1, 'on');

% Convert time vector from seconds to minutes for better readability
time_in_minutes = full_time_out_nb / 60;

% Plot the full deviation over time
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
ax2 = subplot(2, 2, 2);
hold(ax2, 'on');

% Calculate distance from Enceladus, normalized by its SOI radius
distance_from_enceladus_normalized = vecnorm(r_sc_in_flyby_frame_nb, 2, 2)/r_soi_enceladus;
plot(ax2, distance_from_enceladus_normalized, deviation_km, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2, 'DisplayName', 'Position Deviation');

% Mark periapsis
[min_dist, min_idx] = min(distance_from_enceladus_normalized);
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
ax3 = subplot(2, 2, 3);
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
ax4 = subplot(2, 2, 4);
hold(ax4, 'on');

plot(ax4, distance_from_enceladus_normalized, sma_difference_km, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2, 'DisplayName', 'SMA Difference');

% Mark periapsis
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
% createAnimatedFlybyPlot(full_state_out_nb, state_LC_resampled, r_enc_history, full_time_out_nb, pars);

% --- Animated Enceladus-Centric Flyby Plot ---
% createAnimatedEnceladusCentricPlot(r_sc_in_flyby_frame_nb, r_sc_in_flyby_frame_LC, r_saturn_in_flyby_frame_nb, full_time_out_nb, pars);


%% ========================================================================
%  10. ANALYZE ROTATION MATRIX DISCONTINUITIES
%  ========================================================================

% This section identifies discontinuities and saves the corresponding
% rotation matrices into a new array for manual inspection.

fprintf('Analyzing rotation matrix discontinuities...\n');

% --- Calculate the step-to-step difference for both matrix histories ---
num_steps = size(Rm_nb_history, 1);
time_in_minutes = full_time_out_nb / 60; % Ensure time vector is available
discontinuity_nb = zeros(num_steps - 1, 1);
discontinuity_lc = zeros(num_steps - 1, 1);

for i = 2:num_steps
    % Extract consecutive matrices for the N-Body history
    Rm_nb_current = squeeze(Rm_nb_history(i, :, :));
    Rm_nb_previous = squeeze(Rm_nb_history(i-1, :, :));
    
    % Extract consecutive matrices for the Linked Conic history
    Rm_lc_current = squeeze(Rm_LC_history(i, :, :));
    Rm_lc_previous = squeeze(Rm_LC_history(i-1, :, :));

    % Calculate the Frobenius norm of the difference
    discontinuity_nb(i-1) = norm(Rm_nb_current - Rm_nb_previous, 'fro');
    discontinuity_lc(i-1) = norm(Rm_lc_current - Rm_lc_previous, 'fro');
end

% --- Create a new figure for the discontinuity plots ---
figure('Name', 'Rotation Matrix Discontinuity Analysis', 'Color', 'w', 'Position', [200 200 1200 700]);

% Plot for N-Body Rotation Matrix (Rm_nb)
ax_nb_disc = subplot(2, 1, 1);
plot(ax_nb_disc, time_in_minutes(2:end), discontinuity_nb, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.5);
title(ax_nb_disc, 'Discontinuity in N-Body Rotation Matrix (Rm_{nb})');
ylabel(ax_nb_disc, 'Frobenius Norm of Difference');
xlabel(ax_nb_disc, 'Time from Periapsis [minutes]');
grid(ax_nb_disc, 'on');
xlim(ax_nb_disc, [time_in_minutes(2), time_in_minutes(end)]);

% Plot for Linked Conic Rotation Matrix (Rm_LC)
ax_lc_disc = subplot(2, 1, 2);
plot(ax_lc_disc, time_in_minutes(2:end), discontinuity_lc, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.5);
title(ax_lc_disc, 'Discontinuity in Linked Conic Rotation Matrix (Rm_{LC})');
ylabel(ax_lc_disc, 'Frobenius Norm of Difference');
xlabel(ax_lc_disc, 'Time from Periapsis [minutes]');
grid(ax_lc_disc, 'on');
xlim(ax_lc_disc, [time_in_minutes(2), time_in_minutes(end)]);

% Link the x-axes
linkaxes([ax_nb_disc, ax_lc_disc], 'x');


% --- START: AUTOMATIC SPIKE DETECTION ---

% Define a threshold to identify a "big" discontinuity.
discontinuity_threshold = 0.01; 
min_peak_dist_indices = 50; % Minimum distance between peaks in terms of array indices

% Find the peaks (spikes) in the discontinuity data
[pks, locs] = findpeaks(discontinuity_nb, ...
    'MinPeakHeight', discontinuity_threshold, ...
    'MinPeakDistance', min_peak_dist_indices);

% Check if any spikes were found
if ~isempty(locs)
    % The discontinuity at index 'i' corresponds to the jump TO time at 'i+1'
    spike_times_in_minutes = time_in_minutes(locs + 1);
    
    fprintf('\nDetected %d major discontinuity spikes at times [minutes]:\n', length(spike_times_in_minutes));
    fprintf('%.2f\n', spike_times_in_minutes);
    
    % Add vertical lines to the discontinuity plot to mark the detected spikes
    hold(ax_nb_disc, 'on');
    plot(ax_nb_disc, spike_times_in_minutes, pks, 'mv', 'MarkerFaceColor', 'm', 'DisplayName', 'Detected Spikes');
    hold(ax_nb_disc, 'off');
    legend(ax_nb_disc, 'show');

else
    fprintf('\nNo major discontinuities found with threshold > %.2f.\n', discontinuity_threshold);
    spike_times_in_minutes = []; % Ensure variable exists
end

fprintf('Discontinuity analysis complete.\n');


%% ========================================================================
%  11. PLOT ANGULAR MOMENTUM IN THE FLYBY BODY FRAME
%  ========================================================================

% Calculate the specific angular momentum vector (h = r x v) for each time step
% in the Enceladus-centered flyby frame. The cross product is calculated row-wise.
h_flyby_body = cross(r_sc_in_flyby_frame_nb, v_sc_in_flyby_frame_nb, 2);

% Create a new figure for the angular momentum plots
figure('Name', 'Specific Angular Momentum in Flyby Frame', 'Color', 'w', 'Position', [300 300 1000 800]);
sgtitle('Specific Angular Momentum (h = r_{sc} x v_{flyby\_body}) Evolution', 'FontWeight', 'bold');

% --- SUBPLOT 1: Hx component vs. Time ---
ax_hx = subplot(3, 1, 1);
plot(ax_hx, time_in_minutes, h_flyby_body(:,1), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
title(ax_hx, 'X-Component of Angular Momentum (h_x)');
ylabel(ax_hx, 'h_x [km^2/s]');
grid(ax_hx, 'on');
xline(ax_hx, 0, 'r--', 'LineWidth', 1, 'Label', 'Periapsis');

% --- SUBPLOT 2: Hy component vs. Time ---
ax_hy = subplot(3, 1, 2);
plot(ax_hy, time_in_minutes, h_flyby_body(:,2), 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2);
title(ax_hy, 'Y-Component of Angular Momentum (h_y)');
ylabel(ax_hy, 'h_y [km^2/s]');
grid(ax_hy, 'on');
xline(ax_hy, 0, 'r--', 'LineWidth', 1, 'Label', 'Periapsis');

% --- SUBPLOT 3: Hz component vs. Time ---
ax_hz = subplot(3, 1, 3);
plot(ax_hz, time_in_minutes, h_flyby_body(:,3), 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2);
title(ax_hz, 'Z-Component of Angular Momentum (h_z)');
xlabel(ax_hz, 'Time from Periapsis [minutes]');
ylabel(ax_hz, 'h_z [km^2/s]');
grid(ax_hz, 'on');
xline(ax_hz, 0, 'r--', 'LineWidth', 1, 'Label', 'Periapsis');

% Link the x-axes for synchronized zooming and panning
linkaxes([ax_hx, ax_hy, ax_hz], 'x');
xlim(ax_hx, [time_in_minutes(1), time_in_minutes(end)]);


%% ========================================================================
%  12. ANALYZE ROTATION MATRIX BASIS VECTORS
%  ========================================================================

% --- Extract the b2 and b3 basis vectors from the Rm_nb history ---
% Since Rm = [b1'; b2'; b3'], the basis vectors are the rows of the matrix.
b2_history = squeeze(Rm_nb_history(:, 2, :));
b3_history = squeeze(Rm_nb_history(:, 3, :));

% --- Create a new figure to plot the components of b2 and b3 ---
figure('Name', 'Flyby Frame Basis Vector Components', 'Color', 'w', 'Position', [400 400 1200 800]);
sgtitle('Evolution of b_2 and b_3 Basis Vector Components (N-Body)', 'FontWeight', 'bold');

% --- Plot components of b2 ---
ax_b2x = subplot(2, 3, 1); plot(ax_b2x, time_in_minutes, b2_history(:,1), 'b-'); title(ax_b2x, 'b_2 X-Component'); ylabel(ax_b2x, 'b_{2x}'); grid on; xline(ax_b2x, 0, 'r--');
ax_b2y = subplot(2, 3, 2); plot(ax_b2y, time_in_minutes, b2_history(:,2), 'b-'); title(ax_b2y, 'b_2 Y-Component'); ylabel(ax_b2y, 'b_{2y}'); grid on; xline(ax_b2y, 0, 'r--');
ax_b2z = subplot(2, 3, 3); plot(ax_b2z, time_in_minutes, b2_history(:,3), 'b-'); title(ax_b2z, 'b_2 Z-Component'); ylabel(ax_b2z, 'b_{2z}'); grid on; xline(ax_b2z, 0, 'r--');

% --- Plot components of b3 ---
ax_b3x = subplot(2, 3, 4); plot(ax_b3x, time_in_minutes, b3_history(:,1), 'g-'); title(ax_b3x, 'b_3 X-Component'); ylabel(ax_b3x, 'b_{3x}'); grid on; xlabel(ax_b3x, 'Time [min]'); xline(ax_b3x, 0, 'r--');
ax_b3y = subplot(2, 3, 5); plot(ax_b3y, time_in_minutes, b3_history(:,2), 'g-'); title(ax_b3y, 'b_3 Y-Component'); ylabel(ax_b3y, 'b_{3y}'); grid on; xlabel(ax_b3y, 'Time [min]'); xline(ax_b3y, 0, 'r--');
ax_b3z = subplot(2, 3, 6); plot(ax_b3z, time_in_minutes, b3_history(:,3), 'g-'); title(ax_b3z, 'b_3 Z-Component'); ylabel(ax_b3z, 'b_{3z}'); grid on; xlabel(ax_b3z, 'Time [min]'); xline(ax_b3z, 0, 'r--');

% --- Link all axes and add spike lines ---
all_axes_b = [ax_b2x, ax_b2y, ax_b2z, ax_b3x, ax_b3y, ax_b3z];
linkaxes(all_axes_b, 'x');
xlim(ax_b2x, [time_in_minutes(1), time_in_minutes(end)]);
for i = 1:length(all_axes_b)
    for k = 1:length(spike_times_in_minutes)
        xline(all_axes_b(i), spike_times_in_minutes(k), '--m', 'HandleVisibility', 'off');
    end
end


%% ========================================================================
%  13. PLOT r_sc and v_flyby_body over time
%  ========================================================================

% These are the raw vectors in the Saturn-centric J2000 frame that are
% used as inputs to the buildRm function.
r_sc_history = full_state_out_nb(:, 1:3);
% v_flyby_body is already stored in v_enc_history

figure('Name', 'Inertial State Vectors Over Time', 'Color', 'w', 'Position', [500 500 1200 800]);
sgtitle('Evolution of Input Vectors for buildRm (Inertial J2000 Frame)', 'FontWeight', 'bold');

% --- Plot components of r_sc ---
ax_rx = subplot(2, 3, 1); plot(ax_rx, time_in_minutes, r_sc_history(:,1), 'Color', [0.6350 0.0780 0.1840]); title(ax_rx, 'r_{sc} X-Component'); ylabel(ax_rx, 'X [km]'); grid on; xline(ax_rx, 0, 'r--');
ax_ry = subplot(2, 3, 2); plot(ax_ry, time_in_minutes, r_sc_history(:,2), 'Color', [0.6350 0.0780 0.1840]); title(ax_ry, 'r_{sc} Y-Component'); ylabel(ax_ry, 'Y [km]'); grid on; xline(ax_ry, 0, 'r--');
ax_rz = subplot(2, 3, 3); plot(ax_rz, time_in_minutes, r_sc_history(:,3), 'Color', [0.6350 0.0780 0.1840]); title(ax_rz, 'r_{sc} Z-Component'); ylabel(ax_rz, 'Z [km]'); grid on; xline(ax_rz, 0, 'r--');

% --- Plot components of v_flyby_body ---
ax_vx = subplot(2, 3, 4); plot(ax_vx, time_in_minutes, v_enc_history(:,1), 'Color', [0.3010 0.7450 0.9330]); title(ax_vx, 'v_{flyby\_body} X-Component'); ylabel(ax_vx, 'V_x [km/s]'); grid on; xlabel(ax_vx, 'Time [min]'); xline(ax_vx, 0, 'r--');
ax_vy = subplot(2, 3, 5); plot(ax_vy, time_in_minutes, v_enc_history(:,2), 'Color', [0.3010 0.7450 0.9330]); title(ax_vy, 'v_{flyby\_body} Y-Component'); ylabel(ax_vy, 'V_y [km/s]'); grid on; xlabel(ax_vy, 'Time [min]'); xline(ax_vy, 0, 'r--');
ax_vz = subplot(2, 3, 6); plot(ax_vz, time_in_minutes, v_enc_history(:,3), 'Color', [0.3010 0.7450 0.9330]); title(ax_vz, 'v_{flyby\_body} Z-Component'); ylabel(ax_vz, 'V_z [km/s]'); grid on; xlabel(ax_vz, 'Time [min]'); xline(ax_vz, 0, 'r--');

% --- Link all axes and add spike lines ---
all_axes_inputs = [ax_rx, ax_ry, ax_rz, ax_vx, ax_vy, ax_vz];
linkaxes(all_axes_inputs, 'x');
xlim(ax_rx, [time_in_minutes(1), time_in_minutes(end)]);
for i = 1:length(all_axes_inputs)
    for k = 1:length(spike_times_in_minutes)
        xline(all_axes_inputs(i), spike_times_in_minutes(k), '--m', 'HandleVisibility', 'off');
    end
end

%% ========================================================================
%  14. EXTRACT STATE VECTORS AT FIRST DISCONTINUITY
%  ========================================================================

% This section extracts the raw inertial vectors that caused the first major 
% discontinuity spike for detailed analysis.

% --- Check if any spikes were detected in the previous section ---
if ~isempty(spike_times_in_minutes)
    
    % --- Get the time of the very first spike ---
    first_spike_time_sec = spike_times_in_minutes(1) * 60;
    
    % --- Find the index in the master time array closest to the spike time ---
    % We find the minimum absolute difference to get the nearest time step.
    [~, spike_index] = min(abs(full_time_out_nb - first_spike_time_sec));
    
    % --- Extract the state vectors at that specific index ---
    r_sc_at_first_discontinuity = full_state_out_nb(spike_index, 1:3);
    v_flyby_at_first_discontinuity = v_enc_history(spike_index, :);
    
    % --- Display a confirmation message in the Command Window ---
    fprintf('\n--- State at First Discontinuity (t = %.2f min) --- \n', spike_times_in_minutes(1));
    fprintf('Variables created for inspection:\n');
    fprintf('  - r_sc_at_first_discontinuity:  [%.2f, %.2f, %.2f] km\n', r_sc_at_first_discontinuity);
    fprintf('  - v_flyby_at_first_discontinuity: [%.2f, %.2f, %.2f] km/s\n', v_flyby_at_first_discontinuity);
    
else
    fprintf('\nNo discontinuity spikes were detected, so no state vectors were extracted.\n');
end

%% ========================================================================
%  15. PLOT NORM OF CROSS PRODUCT OF INPUT VECTORS FOR buildRm
%  ========================================================================

% Calculate the cross product of the spacecraft's position (inertial, Saturn-centric)
% and the flyby body's velocity (Enceladus's inertial velocity, Saturn-centric).
% This is directly related to the stability of the 'b3' vector in buildRm.
cross_product_input = cross(r_sc_history, v_enc_history, 2); % The '2' ensures row-wise cross product

% Calculate the norm (magnitude) of each resulting vector
norm_cross_product_input = vecnorm(cross_product_input, 2, 2);

figure('Name', 'Norm of Cross Product (r_sc x v_flyby_body)', 'Color', 'w');
ax_norm_cross = axes; % Create axes object for easier adding of lines

hold(ax_norm_cross, 'on');
plot(ax_norm_cross, time_in_minutes, norm_cross_product_input, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2);
title(ax_norm_cross, 'Norm of (r_{sc} x v_{flyby\_body})');
xlabel(ax_norm_cross, 'Time from Periapsis [minutes]');
ylabel(ax_norm_cross, 'Magnitude [km^2/s]');
grid(ax_norm_cross, 'on');
xline(ax_norm_cross, 0, 'r--', 'LineWidth', 1, 'Label', 'Periapsis');

% Add vertical lines for detected discontinuity spikes
if ~isempty(spike_times_in_minutes)
    for k = 1:length(spike_times_in_minutes)
        xline(ax_norm_cross, spike_times_in_minutes(k), '--m', 'HandleVisibility', 'off');
    end
end
hold(ax_norm_cross, 'off');

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