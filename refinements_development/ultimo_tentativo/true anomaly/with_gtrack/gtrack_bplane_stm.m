%% STM B-Plane Targeting with Three-Trajectory Animation
clc;
clear all;
close all;

%% ========================================================================
% 1. SETUP
% ========================================================================
pars.backward_true_anomaly_deg = 100; % Degrees to propagate backward from pericenter
pars.GroundTr.npoints = 30000; % More points for smoother animation
% pars.INPUTS.perturbingBodyNaifIDs = [-2,10];
% pars.INPUTS.perturbingBodyNaifIDs = [-2,10,602, 606]; % J2, Enceladus, Sun, Titan
pars.INPUTS.perturbingBodyNaifIDs = []; % J2, Sun, Enceladus, Titan

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

% Flyby nodes - Partial-COT 1 (O/I) - 1st Flyby
nodein = [4, deg2rad(8.6918), deg2rad(-86.9406)];
nodeout = [4, deg2rad(8.6918), deg2rad(-88.1610)];


spiceParam.frame = 'J2000';
spiceParam.abcorr = 'NONE';
spiceParam.observer = '699';

pars.INPUTS.Flyby.hMapping = 300;

%% ========================================================================
% 2. RUN ENHANCED STM OPTIMIZATION 
% ========================================================================
[NB_in_corrected, B_achieved_final, actual_time_diff, deltaV_magnitude, ...
 iterations_count, converged_flag, B_R_error, B_T_error, backward_duration, ...
 NB_in_state, initial_epoch_nbody, initial_epoch_lc, backward_duration_sec, ...
 mu_TBs, actualBodyNaifIDs, specialPerturbationIDs, flyby_states, ...
 T_linked_conics, propagation_duration, ephem_handles] = ...
 STM_tue_anomaly2(pars, nodein, nodeout);

%% ========================================================================
% 3. SETUP TIME VECTORS FOR PROPAGATION 
% ========================================================================
time_vector_bwd = linspace(0, -propagation_duration*1.2, pars.GroundTr.npoints);
time_vector_fwd = linspace(eps, +propagation_duration*1.2, pars.GroundTr.npoints);

%% ========================================================================
% 4. PROPAGATE ALL THREE TRAJECTORIES - Much cleaner!
% ========================================================================

%% 4a. STM-Corrected N-Body Trajectory
[time_stm_bwd, state_stm_bwd, ~] = propagateNBodyWithSTM(NB_in_corrected(1:3), NB_in_corrected(4:6), ...
    time_vector_bwd, pars.Planet.mu, mu_TBs, ephem_handles.nb_bodies, specialPerturbationIDs, []);

[time_stm_fwd, state_stm_fwd, ~] = propagateNBodyWithSTM(NB_in_corrected(1:3), NB_in_corrected(4:6), ...
    time_vector_fwd, pars.Planet.mu, mu_TBs, ephem_handles.nb_bodies, specialPerturbationIDs, []);

% Create the time-shifted merged trajectory
time_shift = -backward_duration_sec; % This is te_bwd from the STM function

time_stm_bwd_reversed = flip(time_stm_bwd);
state_stm_bwd_reversed = flip(state_stm_bwd, 1);

time_stm_original = [time_stm_bwd_reversed; time_stm_fwd];
state_stm_original = [state_stm_bwd_reversed; state_stm_fwd];

time_stm_merged = time_stm_original - time_shift;

% Interpolate states to match the shifted time grid
state_stm_merged = zeros(size(state_stm_original));
for i = 1:6  % 6 components: 3 position + 3 velocity
    state_stm_merged(:,i) = interp1(time_stm_original, state_stm_original(:,i), ...
                                   time_stm_merged, 'spline', 'extrap');
end

%% 4b. Original N-Body Trajectory (uncorrected) 
[time_orig_bwd, state_orig_bwd, ~] = propagateNBodyWithSTM(NB_in_state(1:3), NB_in_state(4:6), ...
    time_vector_bwd, pars.Planet.mu, mu_TBs, ephem_handles.nb_bodies, specialPerturbationIDs, []);

[time_orig_fwd, state_orig_fwd, ~] = propagateNBodyWithSTM(NB_in_state(1:3), NB_in_state(4:6), ...
     time_vector_fwd, pars.Planet.mu, mu_TBs, ephem_handles.nb_bodies, specialPerturbationIDs, []);

% Create the time-shifted merged trajectory
time_orig_bwd_reversed = flip(time_orig_bwd);
state_orig_bwd_reversed = flip(state_orig_bwd, 1);

time_orig_original = [time_orig_bwd_reversed; time_orig_fwd];
state_orig_original = [state_orig_bwd_reversed; state_orig_fwd];

time_orig_merged = time_orig_original - time_shift;

% Interpolate states to match the shifted time grid
state_orig_merged = zeros(size(state_orig_original));
for i = 1:6  % 6 components: 3 position + 3 velocity
    state_orig_merged(:,i) = interp1(time_orig_original, state_orig_original(:,i), ...
                                    time_orig_merged, 'spline', 'extrap');
end

%% 4c. Linked Conics Trajectory - Using flyby states from STM function
[time_lc_bwd, state_lc_bwd, ~, ~, ~] = propagateKeplerODE2(flyby_states.r0_sc_in', flyby_states.v0_sc_in', ...
    time_vector_bwd, pars.Planet.mu, []);
[time_lc_fwd, state_lc_fwd, ~, ~, ~] = propagateKeplerODE2(flyby_states.r0_sc_out', flyby_states.v0_sc_out', ...
    time_vector_fwd, pars.Planet.mu, []);

% Merge together
time_lc_bwd_reversed = flip(time_lc_bwd);
state_lc_bwd_reversed = flip(state_lc_bwd, 1);

time_lc_merged = [time_lc_bwd_reversed; time_lc_fwd];
state_lc_merged = [state_lc_bwd_reversed; state_lc_fwd];

%% ========================================================================
% 5. GENERATE MOON EPHEMERIS FOR ANIMATION 
% ========================================================================

% Enceladus positions - using the pre-made ephemeris handle
r_enc_history = zeros(length(time_lc_merged), 3);
for i = 1:length(time_lc_merged)
    r_enc = ephem_handles.enceladus_lc(time_lc_merged(i));
    r_enc_history(i, :) = r_enc';
end

% Titan positions - using the pre-made ephemeris handle
r_titan_history = zeros(length(time_lc_merged), 3);
for i = 1:length(time_lc_merged)
    r_titan = ephem_handles.titan_lc(time_lc_merged(i));
    r_titan_history(i, :) = r_titan';
end

%% ========================================================================
% 6. PLOT TRAJECTORIES
% ========================================================================

liveplot_three_trajectories(state_lc_merged, state_stm_merged,...
    state_orig_merged, r_enc_history, r_titan_history, time_lc_merged, pars);

% simple_plot_linked_conics(state_lc_merged, r_enc_history, r_titan_history, pars)

%% ========================================================================
% 7. LINKED CONICS WITH PERICENTRE ORBIT ANALYSIS
% ========================================================================

fprintf('\n=== LINKED CONICS PERICENTRE ORBIT ANALYSIS ===\n');

%% 7a. Calculate Pericentre State Using Flyby Geometry
% Use the same approach as in your reference code
[r_enceladus_at_flyby, v_enceladus_at_flyby] = EphSS_car_spice2(602, pars.INPUTS.epoch0, true, spiceParam);
[vvinfin, r0_sc_in, v0_sc_in, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodein(1), nodein(2), nodein(3), pars.INPUTS.epoch0, pars);
[vvinfout, r0_sc_out, v0_sc_out, vvga] = vinfAlphaCrank_to_VinfCARTClaudio(nodeout(1), nodeout(2), nodeout(3), pars.INPUTS.epoch0, pars);

rp_flyby = pars.INPUTS.Flyby.min_h + pars.Moon.EquRad;
e_fly = 1 + ((rp_flyby*pars.INPUTS.V_inf^2)/pars.Moon.mu);
delta_max = 2*asin(1/e_fly);
pars.delta_max = delta_max;

% Build rotation matrix
b1 = -r0_sc_out./norm(r0_sc_out);
b3 = cross(r0_sc_out, vvga)./norm(cross(r0_sc_out, vvga));
b2 = cross(b3,b1);
Rm = [ b1' b2' b3' ]';

[vvinfouBM, delta] = infVelBeforeDefect(vvinfin, vvinfout, delta_max, vvga);
vvinfin_bf = [Rm*vvinfin']';
vvinfouBM_bf = [Rm*vvinfouBM']';

% Calculate orbital parameters
Energy = 0.5*norm(vvinfin_bf)^2;
sma = -pars.Moon.mu/(2*Energy);
ecc = 1/(sin(delta/2));
rp = sma*(1 - ecc);
hhat = cross(vvinfin_bf, vvinfouBM_bf)./norm(cross(vvinfin_bf, vvinfouBM_bf));
vp = sqrt(norm(vvinfin_bf)^2 + 2*pars.Moon.mu/rp);

% Calculate pericentre state in flyby frame
rrp_bf = rp.*(vvinfin_bf - vvinfouBM_bf)./norm(vvinfin_bf - vvinfouBM_bf);
vvp_bf = vp.*cross(hhat, rrp_bf./rp);

% Transform to Saturn-centric frame
Rm_inv = Rm';
rrp_saturn_centric_lc = (Rm_inv * rrp_bf')' + r_enceladus_at_flyby;
vvp_saturn_centric_lc = (Rm_inv * vvp_bf')' + v_enceladus_at_flyby;

rrp_saturn_centric_lc = rrp_saturn_centric_lc(:);
vvp_saturn_centric_lc = vvp_saturn_centric_lc(:);

fprintf('  Pericentre altitude: %.3f km\n', rp - pars.Moon.EquRad);
fprintf('  Pericentre velocity: %.3f km/s\n', vp);

%% 7b. Propagate Linked Conics Orbit (2-Body Dynamics)
% Use 2-body propagation by setting empty perturbation lists
mu_central_body = pars.Planet.mu;
empty_mu_TBs = [];
empty_ephem_handle = @(t) [];  % Empty ephemeris handle
empty_special_IDs = [];


% Forward propagation
[time_out_fwd_lc, state_out_fwd_lc] = propagateNBodyODE2(rrp_saturn_centric_lc, vvp_saturn_centric_lc, ...
    time_vector_fwd, mu_central_body, empty_mu_TBs, empty_ephem_handle, empty_special_IDs);

% Backward propagation  
[time_out_bwd_lc, state_out_bwd_lc] = propagateNBodyODE2(rrp_saturn_centric_lc, vvp_saturn_centric_lc, ...
    time_vector_bwd, mu_central_body, empty_mu_TBs, empty_ephem_handle, empty_special_IDs);

% Merge trajectories
time_out_bwd_lc = flipud(time_out_bwd_lc);
state_out_bwd_lc = flipud(state_out_bwd_lc);
full_time_out_lc = [time_out_bwd_lc; time_out_fwd_lc];
full_state_out_lc = [state_out_bwd_lc; state_out_fwd_lc];
num_steps_lc = length(full_time_out_lc);

%% 7c. Generate Enceladus Ephemeris for Linked Conics Analysis
r_enc_history_lc = zeros(num_steps_lc, 3);
v_enc_history_lc = zeros(num_steps_lc, 3);

fprintf('  Generating Enceladus ephemeris...\n');

for i = 1:num_steps_lc
    current_mjd = pars.INPUTS.epoch0 + full_time_out_lc(i) / 86400;
    [r_enc, v_enc] = EphSS_car_spice2(602, current_mjd, true, spiceParam);
    r_enc_history_lc(i, :) = r_enc;
    v_enc_history_lc(i, :) = v_enc;
end

%% 7d. Transform to Enceladus Body-Fixed Frame and Calculate Ground Track
lats_lc = zeros(num_steps_lc, 1);
lon_rad_temp_lc = zeros(num_steps_lc, 1);
alts_lc = zeros(num_steps_lc, 1);
r_in_flyby_frame_lc = zeros(num_steps_lc, 3);

fprintf('  Computing ground track coordinates...\n');

for i = 1:num_steps_lc
    current_mjd = pars.INPUTS.epoch0 + full_time_out_lc(i) / 86400;
    sc_pos_icrf = full_state_out_lc(i, 1:3);
    r_enc = r_enc_history_lc(i, :);
    v_enc = v_enc_history_lc(i, :);
    
    % --- TRANSFORMATION USING time-varying Rm ---
    sc_pos_relative_icrf = sc_pos_icrf - r_enc;
    Rm_current = buildRm(r_enc, v_enc);
    r_in_flyby_frame_lc(i, :) = (Rm_current * sc_pos_relative_icrf')';
    r_in_flyby_frame_unit = r_in_flyby_frame_lc(i, :) / norm(r_in_flyby_frame_lc(i, :));
    
    lon_rad_temp_lc(i) = atan2(r_in_flyby_frame_unit(2), r_in_flyby_frame_unit(1));
    lats_lc(i) = asin(r_in_flyby_frame_unit(3));
    
    % Altitude calculation
    alts_lc(i) = norm(sc_pos_relative_icrf) - pars.Moon.EquRad;
end

% Convert to east-longitude format
long_lc = wrapTo2Pi(lon_rad_temp_lc);
longs_lc = wrapTo360(long_lc);

% Convert to degrees
lats_deg_lc = rad2deg(lats_lc);
longs_deg_lc = rad2deg(longs_lc);

% Find periapsis index
[~, idx_periapsis_lc] = min(alts_lc);

% Create data structure for plotting
mapping_altitude_km = pars.INPUTS.Flyby.hMapping;
valid_indices_lc = alts_lc <= mapping_altitude_km;

LinkedConics_Mapping_Data.lats = lats_lc(valid_indices_lc);
LinkedConics_Mapping_Data.longs = longs_lc(valid_indices_lc);
LinkedConics_Mapping_Data.rp_lat = lats_deg_lc(idx_periapsis_lc);
LinkedConics_Mapping_Data.rp_long = longs_deg_lc(idx_periapsis_lc);

fprintf('  Pericentre ground track: Lat = %.3f°, Lon = %.3f°\n', ...
    LinkedConics_Mapping_Data.rp_lat, LinkedConics_Mapping_Data.rp_long);

%% 7e. Plot Linked Conics Analysis Results
fprintf('  Creating plots...\n');

% --- Create comprehensive figure ---
figure('Name', 'Linked Conics Pericentre Orbit Analysis', 'Color', 'w', 'Position', [150 150 1600 900]);

% --- SUBPLOT 1: 3D Trajectory Comparison ---
ax1 = subplot(2, 2, 1);
hold(ax1, 'on');

% Plot trajectories
plot3(ax1, full_state_out_lc(:,1), full_state_out_lc(:,2), full_state_out_lc(:,3), ...
    'g-', 'LineWidth', 1, 'DisplayName', 'Linked Conics (Pericentre)');
plot3(ax1, r_enc_history_lc(:,1), r_enc_history_lc(:,2), r_enc_history_lc(:,3), ...
    'k-', 'LineWidth', 1, 'DisplayName', 'Enceladus');

% Mark key points
plot3(ax1, rrp_saturn_centric_lc(1), rrp_saturn_centric_lc(2), rrp_saturn_centric_lc(3), ...
    'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red', 'DisplayName', 'Pericentre');

title(ax1, '3D Trajectory - Linked Conics from Pericentre');
xlabel(ax1, 'X [km]'); ylabel(ax1, 'Y [km]'); zlabel(ax1, 'Z [km]');
legend(ax1, 'show', 'Location', 'best');
grid(ax1, 'on'); axis(ax1, 'equal'); view(ax1, 3);

% --- SUBPLOT 2: Ground Track ---
ax2 = subplot(2, 2, 2);
hold(ax2, 'on');

% Plot texture map
plotTextureLatLong(pars.INPUTS.idMoon, pars.INPUTS.idCentral, 1);

% Plot ground track
if ~isempty(LinkedConics_Mapping_Data.lats)
    Plot_Flyby_GTClaudio(LinkedConics_Mapping_Data, [0.2, 0.8, 0.2], '-', 3);
end

title(ax2, sprintf('Ground Track (Altitude ≤ %d km)', mapping_altitude_km));
grid(ax2, 'on');

% --- SUBPLOT 3: Altitude vs Time ---
ax3 = subplot(2, 2, 3);
time_minutes_lc = full_time_out_lc / 60; % CONVERT TO MINUTES
plot(ax3, time_minutes_lc, alts_lc, 'g-', 'LineWidth', 1);
hold(ax3, 'on');

% Mark periapsis
plot(ax3, time_minutes_lc(idx_periapsis_lc), alts_lc(idx_periapsis_lc), ...
    'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red');

% Add mapping altitude line
yline(ax3, mapping_altitude_km, 'r--', 'LineWidth', 1, 'DisplayName', 'Mapping Altitude');

title(ax3, 'Altitude Profile');
xlabel(ax3, 'Time from Pericentre [minutes]'); % UPDATE LABEL
ylabel(ax3, 'Altitude [km]');
grid(ax3, 'on');
legend(ax3, 'show');
xlim(ax3, [-10, 10]); % ZOOM AXES

% --- SUBPLOT 4: Trajectory in Enceladus Frame ---
ax4 = subplot(2, 2, 4);
hold(ax4, 'on');

% Plot Enceladus as sphere
[x_sphere, y_sphere, z_sphere] = sphere(50);
radius = pars.Moon.EquRad;
surf(ax4, x_sphere*radius, y_sphere*radius, z_sphere*radius, ...
     'FaceColor', [0.8, 0.9, 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.7, 'DisplayName', 'Enceladus');

% Plot trajectory in Enceladus frame
plot3(ax4, r_in_flyby_frame_lc(:,1), r_in_flyby_frame_lc(:,2), r_in_flyby_frame_lc(:,3), ...
    'g-', 'LineWidth', 1, 'DisplayName', 'Linked Conics');

% Mark pericentre
plot3(ax4, r_in_flyby_frame_lc(idx_periapsis_lc,1), ...
          r_in_flyby_frame_lc(idx_periapsis_lc,2), ...
          r_in_flyby_frame_lc(idx_periapsis_lc,3), ...
    'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red', 'DisplayName', 'Pericentre');

title(ax4, 'Trajectory in Enceladus-Centered Frame');
xlabel(ax4, 'X [km]'); ylabel(ax4, 'Y [km]'); zlabel(ax4, 'Z [km]');
legend(ax4, 'show');
grid(ax4, 'on'); axis(ax4, 'equal'); view(ax4, 3);

