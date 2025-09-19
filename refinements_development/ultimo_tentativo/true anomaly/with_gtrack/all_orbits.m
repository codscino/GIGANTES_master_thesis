%% STM B-Plane Targeting with Three-Trajectory Animation and Comparison
clc;
clear all;
close all;

%% ========================================================================
% 1. SETUP
% ========================================================================
pars.backward_true_anomaly_deg = 80; % Degrees to propagate backward from pericenter
pars.GroundTr.npoints = 30000; % More points for smoother animation
% pars.INPUTS.perturbingBodyNaifIDs = [-2,10];
% pars.INPUTS.perturbingBodyNaifIDs = [-2,10,602, 606]; % J2, Enceladus, Sun, Titan
pars.INPUTS.perturbingBodyNaifIDs = [-2,10,602,606]; % J2, Sun, Enceladus, Titan

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

%% 4b. Original N-Body Trajectory (uncorrected) - No recalculation needed!
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
% 7. THREE-TRAJECTORY COMPARISON ANALYSIS WITH PERICENTRE ORBIT
% ========================================================================

fprintf('\n=== THREE-TRAJECTORY COMPARISON ANALYSIS ===\n');

%% 7a. Calculate Pericentre State Using Flyby Geometry (for Linked Conics)
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
mu_central_body = pars.Planet.mu;
empty_mu_TBs = [];
empty_ephem_handle = @(t) [];
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
num_steps = length(full_time_out_lc);

%% 7c. Generate Enceladus Ephemeris for All Trajectories
r_enc_history_analysis = zeros(num_steps, 3);
v_enc_history_analysis = zeros(num_steps, 3);

fprintf('  Generating Enceladus ephemeris...\n');

for i = 1:num_steps
    current_mjd = pars.INPUTS.epoch0 + full_time_out_lc(i) / 86400;
    [r_enc, v_enc] = EphSS_car_spice2(602, current_mjd, true, spiceParam);
    r_enc_history_analysis(i, :) = r_enc;
    v_enc_history_analysis(i, :) = v_enc;
end

%% 7d. Process All Three Trajectories
fprintf('  Processing all three trajectories...\n');

% Define the three trajectories to process
trajectories = struct();
trajectories.LinkedConics.name = 'Linked Conics';
trajectories.LinkedConics.states = full_state_out_lc;
trajectories.LinkedConics.color = [0.2, 0.8, 0.2]; % Green
trajectories.LinkedConics.linestyle = '-';
trajectories.LinkedConics.linewidth = 3; % Thicker for ground track

trajectories.STMCorrected.name = 'STM-Corrected N-Body';
trajectories.STMCorrected.states = state_stm_merged;
trajectories.STMCorrected.color = [0.8, 0.2, 0.2]; % Red
trajectories.STMCorrected.linestyle = '--'; % Dashed
trajectories.STMCorrected.linewidth = 2;

trajectories.OriginalNBody.name = 'Original N-Body';
trajectories.OriginalNBody.states = state_orig_merged;
trajectories.OriginalNBody.color = [0.2, 0.2, 0.8]; % Blue
trajectories.OriginalNBody.linestyle = ':';
trajectories.OriginalNBody.linewidth = 2;

trajectory_names = fieldnames(trajectories);
mapping_altitude_km = pars.INPUTS.Flyby.hMapping;

% Process each trajectory
for t = 1:length(trajectory_names)
    traj_name = trajectory_names{t};
    current_traj = trajectories.(traj_name);
    
    fprintf('    Processing %s trajectory...\n', current_traj.name);
    
    % Initialize arrays
    lats = zeros(num_steps, 1);
    lon_rad_temp = zeros(num_steps, 1);
    alts = zeros(num_steps, 1);
    r_in_flyby_frame = zeros(num_steps, 3);
    
    % Transform to Enceladus body-fixed frame and calculate ground track
    for i = 1:num_steps
        sc_pos_icrf = current_traj.states(i, 1:3);
        r_enc = r_enc_history_analysis(i, :);
        v_enc = v_enc_history_analysis(i, :);
        
        % Transformation using time-varying Rm
        sc_pos_relative_icrf = sc_pos_icrf - r_enc;
        Rm_current = buildRm(r_enc, v_enc);
        r_in_flyby_frame(i, :) = (Rm_current * sc_pos_relative_icrf')';
        r_in_flyby_frame_unit = r_in_flyby_frame(i, :) / norm(r_in_flyby_frame(i, :));
        
        lon_rad_temp(i) = atan2(r_in_flyby_frame_unit(2), r_in_flyby_frame_unit(1));
        lats(i) = asin(r_in_flyby_frame_unit(3));
        
        % Altitude calculation
        alts(i) = norm(sc_pos_relative_icrf) - pars.Moon.EquRad;
    end
    
    % Convert to east-longitude format
    long = wrapTo2Pi(lon_rad_temp);
    longs = wrapTo360(long);
    
    % Convert to degrees
    lats_deg = rad2deg(lats);
    longs_deg = rad2deg(longs);
    
    % Find periapsis index
    [~, idx_periapsis] = min(alts);
    
    % Store data for this trajectory
    trajectories.(traj_name).lats = lats;
    trajectories.(traj_name).longs = longs;
    trajectories.(traj_name).lats_deg = lats_deg;
    trajectories.(traj_name).longs_deg = longs_deg;
    trajectories.(traj_name).alts = alts;
    trajectories.(traj_name).r_in_flyby_frame = r_in_flyby_frame;
    trajectories.(traj_name).idx_periapsis = idx_periapsis;
    trajectories.(traj_name).rp_lat = lats_deg(idx_periapsis);
    trajectories.(traj_name).rp_long = longs_deg(idx_periapsis);
    
    % Create mapping data (altitude below mapping threshold)
    valid_indices = alts <= mapping_altitude_km;
    trajectories.(traj_name).mapping_data.lats = lats(valid_indices);
    trajectories.(traj_name).mapping_data.longs = longs(valid_indices);
    trajectories.(traj_name).mapping_data.rp_lat = lats_deg(idx_periapsis);
    trajectories.(traj_name).mapping_data.rp_long = longs_deg(idx_periapsis);
    
    fprintf('      Pericentre: Lat = %.3f°, Lon = %.3f°, Alt = %.3f km\n', ...
        trajectories.(traj_name).rp_lat, trajectories.(traj_name).rp_long, alts(idx_periapsis));
end

%% 7e. Plot Three-Trajectory Comparison Results
fprintf('  Creating comparison plots...\n');

% Create comprehensive figure
figure('Name', 'Three-Trajectory Comparison Analysis', 'Color', 'w', 'Position', [150 150 1600 900]);

% --- SUBPLOT 1: 3D Trajectory Comparison ---
ax1 = subplot(2, 2, 1);
hold(ax1, 'on');

% Plot Enceladus orbit
plot3(ax1, r_enc_history_analysis(:,1), r_enc_history_analysis(:,2), r_enc_history_analysis(:,3), ...
    'k-', 'LineWidth', 1, 'DisplayName', 'Enceladus');

% Plot all three trajectories
for t = 1:length(trajectory_names)
    traj_name = trajectory_names{t};
    current_traj = trajectories.(traj_name);
    
    plot3(ax1, current_traj.states(:,1), current_traj.states(:,2), current_traj.states(:,3), ...
        'Color', current_traj.color, 'LineStyle', current_traj.linestyle, 'LineWidth', 1.5, ...
        'DisplayName', current_traj.name);
    
    % Mark pericentre
    idx_p = current_traj.idx_periapsis;
    plot3(ax1, current_traj.states(idx_p,1), current_traj.states(idx_p,2), current_traj.states(idx_p,3), ...
        'o', 'MarkerSize', 6, 'MarkerFaceColor', current_traj.color, 'MarkerEdgeColor', 'black', 'HandleVisibility', 'off');
end

title(ax1, '3D Trajectory Comparison');
xlabel(ax1, 'X [km]'); ylabel(ax1, 'Y [km]'); zlabel(ax1, 'Z [km]');
legend(ax1, 'show', 'Location', 'best');
grid(ax1, 'on'); axis(ax1, 'equal'); view(ax1, 3);

% --- SUBPLOT 2: Ground Track Comparison ---
ax2 = subplot(2, 2, 2);
hold(ax2, 'on');

% Plot texture map
plotTextureLatLong(pars.INPUTS.idMoon, pars.INPUTS.idCentral, 1);

% Plot ground tracks for all trajectories
for t = 1:length(trajectory_names)
    traj_name = trajectory_names{t};
    current_traj = trajectories.(traj_name);
    
    if ~isempty(current_traj.mapping_data.lats)
        Plot_Flyby_GTClaudio_Enhanced(current_traj.mapping_data, current_traj.color, current_traj.name, current_traj.linewidth, current_traj.linestyle);
    end
end

title(ax2, sprintf('Ground Track Comparison (Altitude ≤ %d km)', mapping_altitude_km));
legend(ax2, 'show', 'Location', 'best');
grid(ax2, 'on');

% --- SUBPLOT 3: Altitude vs Time Comparison ---
ax3 = subplot(2, 2, 3);
hold(ax3, 'on');

time_minutes = full_time_out_lc / 60;

for t = 1:length(trajectory_names)
    traj_name = trajectory_names{t};
    current_traj = trajectories.(traj_name);
    
    plot(ax3, time_minutes, current_traj.alts, 'Color', current_traj.color, ...
        'LineStyle', current_traj.linestyle, 'LineWidth', 1.5, 'DisplayName', current_traj.name);
    
    % Mark periapsis
    idx_p = current_traj.idx_periapsis;
    plot(ax3, time_minutes(idx_p), current_traj.alts(idx_p), 'o', 'MarkerSize', 6, ...
        'MarkerFaceColor', current_traj.color, 'MarkerEdgeColor', 'black', 'HandleVisibility', 'off');
end

% Add mapping altitude line
yline(ax3, mapping_altitude_km, 'r--', 'LineWidth', 1, 'DisplayName', 'Mapping Altitude');

title(ax3, 'Altitude Profile Comparison');
xlabel(ax3, 'Time from Pericentre [minutes]');
ylabel(ax3, 'Altitude [km]');
grid(ax3, 'on');
legend(ax3, 'show');
xlim(ax3, [-10, 10]);

% --- SUBPLOT 4: Trajectory in Enceladus Frame Comparison ---
ax4 = subplot(2, 2, 4);
hold(ax4, 'on');

% Plot Enceladus as sphere
[x_sphere, y_sphere, z_sphere] = sphere(50);
radius = pars.Moon.EquRad;
surf(ax4, x_sphere*radius, y_sphere*radius, z_sphere*radius, ...
     'FaceColor', [0.8, 0.9, 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.7, 'DisplayName', 'Enceladus');

% Plot trajectories in Enceladus frame
for t = 1:length(trajectory_names)
    traj_name = trajectory_names{t};
    current_traj = trajectories.(traj_name);
    
    plot3(ax4, current_traj.r_in_flyby_frame(:,1), current_traj.r_in_flyby_frame(:,2), ...
        current_traj.r_in_flyby_frame(:,3), 'Color', current_traj.color, ...
        'LineStyle', current_traj.linestyle, 'LineWidth', 1.5, 'DisplayName', current_traj.name);
    
    % Mark pericentre
    idx_p = current_traj.idx_periapsis;
    plot3(ax4, current_traj.r_in_flyby_frame(idx_p,1), current_traj.r_in_flyby_frame(idx_p,2), ...
        current_traj.r_in_flyby_frame(idx_p,3), 'o', 'MarkerSize', 6, ...
        'MarkerFaceColor', current_traj.color, 'MarkerEdgeColor', 'black', 'HandleVisibility', 'off');
end

title(ax4, 'Trajectory Comparison in Enceladus-Centered Frame');
xlabel(ax4, 'X [km]'); ylabel(ax4, 'Y [km]'); zlabel(ax4, 'Z [km]');
legend(ax4, 'show');
grid(ax4, 'on'); axis(ax4, 'equal'); view(ax4, 3);

%% 7f. Print Summary Statistics
fprintf('\n=== TRAJECTORY COMPARISON SUMMARY ===\n');
for t = 1:length(trajectory_names)
    traj_name = trajectory_names{t};
    current_traj = trajectories.(traj_name);
    fprintf('%s:\n', current_traj.name);
    fprintf('  Pericentre Altitude: %.3f km\n', current_traj.alts(current_traj.idx_periapsis));
    fprintf('  Pericentre Lat/Lon: %.3f°, %.3f°\n', current_traj.rp_lat, current_traj.rp_long);
    fprintf('  Points below %d km: %d\n', mapping_altitude_km, sum(current_traj.alts <= mapping_altitude_km));
end

%% ========================================================================
% HELPER FUNCTIONS
% ========================================================================

function Plot_Flyby_GTClaudio_Enhanced(Flyby_Data, color, trajectory_name, linewidth, linestyle)
% Enhanced version of Plot_Flyby_GT_claudio with trajectory name, linewidth, and linestyle
    lats = Flyby_Data.lats;
    longs = Flyby_Data.longs;
    rp_lat = Flyby_Data.rp_lat;
    rp_long = Flyby_Data.rp_long;
    
    % Convert to degrees for plotting
    lats_deg = rad2deg(lats);
    longs_deg = rad2deg(longs);
    
    % Find if there is a jump in longitude (passing from near 0 to 360)
    stop_index = find(diff(longs) > 0, 1);
    
    % Plot the ground track
    if ~isempty(stop_index) && ~isequal(stop_index, 1)
        % Split the vector at the negative difference indices
        Long_1 = longs(1:stop_index);
        Lat_1 = lats(1:stop_index);
        Long_2 = longs(stop_index+1:end);
        Lat_2 = lats(stop_index+1:end);
        
        plot(rad2deg(Long_1), rad2deg(Lat_1), 'Color', color, 'LineStyle', linestyle, ...
            'LineWidth', linewidth, 'DisplayName', trajectory_name);
        plot(rad2deg(Long_2), rad2deg(Lat_2), 'Color', color, 'LineStyle', linestyle, ...
            'LineWidth', linewidth, 'HandleVisibility', 'off');
    else
        plot(longs_deg, lats_deg, 'Color', color, 'LineStyle', linestyle, ...
            'LineWidth', linewidth, 'DisplayName', trajectory_name);
    end
    
    % Plot pericentre point
    plot(rp_long, rp_lat, 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', color, ...
        'MarkerSize', 6, 'HandleVisibility', 'off');
end