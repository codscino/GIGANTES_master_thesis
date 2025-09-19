%% STM B-Plane Targeting with Three-Trajectory Animation
clc;
clear all;
close all;

%% ========================================================================
% 1. SETUP
% ========================================================================
pars.backward_true_anomaly_deg = 100; % Degrees to propagate backward from pericenter
pars.GroundTr.npoints = 30000; % More points for smoother animation
pars.INPUTS.perturbingBodyNaifIDs = [];
% pars.INPUTS.perturbingBodyNaifIDs = [-2,10,602, 606]; % J2, Enceladus, Sun
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
% pars.EncPlotSize = 40; % For plotting Enceladus

% Flyby nodes - Partial-COT 1 (O/I) - 1st Flyby
nodein = [4, deg2rad(8.6918), deg2rad(-86.9406)];
nodeout = [4, deg2rad(8.6918), deg2rad(-88.1610)];

%% ========================================================================
% 2. RUN STM OPTIMIZATION
% ========================================================================
[NB_in_corrected, B_achieved_final, actual_time_diff, deltaV_magnitude, ...
 iterations_count, converged_flag, B_R_error, B_T_error, backward_duration] = ...
 STM_tue_anomaly(pars, nodein, nodeout);


%% ========================================================================
% 3. CALCULATE LINKED CONICS PERIOD AND PROPAGATION TIME
% ========================================================================
% Get linked conics states at flyby epoch
[vvinfin, r0_sc_in, v0_sc_in, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodein(1), nodein(2), nodein(3), pars.INPUTS.epoch0, pars);
[vvinfout, r0_sc_out, v0_sc_out, vvga] = vinfAlphaCrank_to_VinfCARTClaudio(nodeout(1), nodeout(2), nodeout(3), pars.INPUTS.epoch0, pars);

% Calculate linked conics period using the flyby state
lc_state_flyby = [r0_sc_in, v0_sc_in];
lc_kep = car2kep(lc_state_flyby, pars.Planet.mu);
T_linked_conics = 2 * pi * sqrt(lc_kep(1)^3 / pars.Planet.mu); % Period in seconds
T_linked_conics_hours = T_linked_conics / 3600;

% Propagation duration: T_LC / 2.1 (both backward and forward)
propagation_duration = T_linked_conics/2;
propagation_duration_hours = propagation_duration / 3600;

fprintf('Linked Conics Period: %.2f hours\n', T_linked_conics_hours);
fprintf('Propagation Duration (T_LC/2.1): %.2f hours each direction\n', propagation_duration_hours);

%% ========================================================================
% 4. PROPAGATE ALL THREE TRAJECTORIES
% ========================================================================

% Setup for propagation
actualBodyNaifIDs = pars.INPUTS.perturbingBodyNaifIDs(pars.INPUTS.perturbingBodyNaifIDs >= 0);
mu_TBs = zeros(1, length(actualBodyNaifIDs));
for i = 1:length(actualBodyNaifIDs)
    id = actualBodyNaifIDs(i);
    switch id
        case 10,  mu_TBs(i) = getAstroConstants('Sun', 'Mu');
        case 5,   mu_TBs(i) = getAstroConstants('Jupiter', 'Mu');
        case 601, [~, mu, ~, ~] = satMoonsConstants(0); mu_TBs(i) = mu;
        case 602, [~, mu, ~, ~] = satMoonsConstants(1); mu_TBs(i) = mu;
        case 603, [~, mu, ~, ~] = satMoonsConstants(2); mu_TBs(i) = mu;
        case 604, [~, mu, ~, ~] = satMoonsConstants(3); mu_TBs(i) = mu;
        case 605, [~, mu, ~, ~] = satMoonsConstants(4); mu_TBs(i) = mu;
        case 606, [~, mu, ~, ~] = satMoonsConstants(5); mu_TBs(i) = mu;
        case 607, mu_TBs(i) = 0.374;
        case 608, mu_TBs(i) = 120.4;
    end
end

specialPerturbationIDs = pars.INPUTS.perturbingBodyNaifIDs(pars.INPUTS.perturbingBodyNaifIDs < 0);


% N-Body trajectories: Start from corrected epoch (earlier than flyby)
initial_epoch_nbody = pars.INPUTS.epoch0 - backward_duration/24; % MJD2000
fprintf('N-Body initial epoch: %.6f MJD2000 (%.2f hours before flyby)\n', initial_epoch_nbody, backward_duration);

% Linked Conics: Reference from flyby epoch
initial_epoch_lc = pars.INPUTS.epoch0; % MJD2000
fprintf('Linked Conics reference epoch: %.6f MJD2000 (flyby epoch)\n', initial_epoch_lc);

% Time vectors for propagation
time_vector_bwd = linspace(0, -propagation_duration, pars.GroundTr.npoints);
time_vector_fwd = linspace(eps, +propagation_duration, pars.GroundTr.npoints);

% Ephemeris functions for N-body (relative to corrected epoch)
spiceParam.frame = 'J2000';
spiceParam.abcorr = 'NONE';
spiceParam.observer = '699';

ephem_handle_nb = @(t_sec) get_body_positions_wrapper(t_sec, initial_epoch_nbody, actualBodyNaifIDs, spiceParam);
ephem_enceladus_nb = @(t) get_body_positions_wrapper(t, initial_epoch_nbody, pars.INPUTS.NAIFMoon, spiceParam);

% Ephemeris functions for Linked Conics (relative to flyby epoch)
ephem_enceladus_lc = @(t) get_body_positions_wrapper(t, initial_epoch_lc, pars.INPUTS.NAIFMoon, spiceParam);

%% 4.0 -theta state
% Get original uncorrected initial state by propagating LC backward
initial_kep_lc = car2kep([r0_sc_in, v0_sc_in], pars.Planet.mu);
initial_true_anomaly_deg = rad2deg(initial_kep_lc(6));

% Create event function for backward propagation
event_func = @(t, x) trueAnomalyBackwardEvent(t, x, pars.Planet.mu, initial_true_anomaly_deg, pars.backward_true_anomaly_deg);

% Propagate backward to get original uncorrected state
[~, ~, te_bwd, ye_bwd, ~] = propagateKeplerODE2(r0_sc_in', v0_sc_in', time_vector_bwd, pars.Planet.mu, event_func);

NB_in_state = ye_bwd(1, :)';
NB_in_time = te_bwd; %seconds

%% 4a. STM-Corrected N-Body Trajectory
[time_stm_bwd, state_stm_bwd, ~] = propagateNBodyWithSTM(NB_in_corrected(1:3), NB_in_corrected(4:6), ...
    time_vector_bwd, pars.Planet.mu, mu_TBs, ephem_handle_nb, specialPerturbationIDs, []);

[time_stm_fwd, state_stm_fwd, ~] = propagateNBodyWithSTM(NB_in_corrected(1:3), NB_in_corrected(4:6), ...
    time_vector_fwd, pars.Planet.mu, mu_TBs, ephem_handle_nb, specialPerturbationIDs, []);

% Create the original (unshifted) merged time vector
time_stm_bwd_reversed = flip(time_stm_bwd);
state_stm_bwd_reversed = flip(state_stm_bwd, 1);

time_stm_original = [time_stm_bwd_reversed; time_stm_fwd];
state_stm_original = [state_stm_bwd_reversed; state_stm_fwd];

% Create the desired shifted time vector
time_stm_merged = time_stm_original - te_bwd;

% Interpolate states to match the shifted time grid
% We need to interpolate from the original time grid to the shifted time grid
state_stm_merged = zeros(size(state_stm_original));
for i = 1:6  % 6 components: 3 position + 3 velocity
    state_stm_merged(:,i) = interp1(time_stm_original, state_stm_original(:,i), ...
                                   time_stm_merged, 'spline', 'extrap');
end

%% 4b. Original N-Body Trajectory (uncorrected) - Same approach
[time_orig_bwd, state_orig_bwd, ~] = propagateNBodyWithSTM(NB_in_state(1:3), NB_in_state(4:6), ...
    time_vector_bwd, pars.Planet.mu, mu_TBs, ephem_handle_nb, specialPerturbationIDs, []);

[time_orig_fwd, state_orig_fwd, ~] = propagateNBodyWithSTM(NB_in_state(1:3), NB_in_state(4:6), ...
     time_vector_fwd, pars.Planet.mu, mu_TBs, ephem_handle_nb, specialPerturbationIDs, []);

% Create the original (unshifted) merged time vector
time_orig_bwd_reversed = flip(time_orig_bwd);
state_orig_bwd_reversed = flip(state_orig_bwd, 1);

time_orig_original = [time_orig_bwd_reversed; time_orig_fwd];
state_orig_original = [state_orig_bwd_reversed; state_orig_fwd];

% Create the desired shifted time vector
time_orig_merged = time_orig_original - te_bwd;

% Interpolate states to match the shifted time grid
state_orig_merged = zeros(size(state_orig_original));
for i = 1:6  % 6 components: 3 position + 3 velocity
    state_orig_merged(:,i) = interp1(time_orig_original, state_orig_original(:,i), ...
                                    time_orig_merged, 'spline', 'extrap');
end

%% 4c. Linked Conics Trajectory - No change needed since no shift applied
[time_lc_bwd, state_lc_bwd, ~, ~, ~] = propagateKeplerODE2(r0_sc_in', v0_sc_in',time_vector_bwd , pars.Planet.mu, []);
[time_lc_fwd, state_lc_fwd, ~, ~, ~] = propagateKeplerODE2(r0_sc_out', v0_sc_out',time_vector_fwd , pars.Planet.mu, []);

% merge together:
time_lc_bwd_reversed = flip(time_lc_bwd);
state_lc_bwd_reversed = flip(state_lc_bwd, 1);

time_lc_merged = [time_lc_bwd_reversed; time_lc_fwd];
state_lc_merged = [state_lc_bwd_reversed; state_lc_fwd];

%% ========================================================================
% 5. GENERATE MOON EPHEMERIS FOR ANIMATION
% ========================================================================

% Enceladus positions (relative to N-body epoch for proper synchronization)
r_enc_history = zeros(length(time_lc_merged), 3);
for i = 1:length(time_lc_merged)
    r_enc = ephem_enceladus_lc(time_lc_merged(i));
    r_enc_history(i, :) = r_enc';
end

% Titan positions (relative to N-body epoch)
r_titan_history = zeros(length(time_lc_merged), 3);
ephem_titan_nb = @(t) get_body_positions_wrapper(t, initial_epoch_lc, 606, spiceParam);
for i = 1:length(time_lc_merged)
    r_titan = ephem_titan_nb(time_lc_merged(i));
    r_titan_history(i, :) = r_titan';
end


%% ========================================================================
% 6. plot trajectories
% ========================================================================

liveplot_three_trajectories(state_lc_merged, state_stm_merged,...
    state_orig_merged, r_enc_history, r_titan_history, time_lc_merged, pars);

% plot_static_three_trajectories(state_lc_merged, state_stm_merged, state_orig_merged, ...
%                               r_enc_history, r_titan_history, time_lc_merged, backward_duration, pars);



%% ========================================================================
% HELPER FUNCTIONS
% ========================================================================

function [value, isterminal, direction] = trueAnomalyBackwardEvent(t, x, mu, initial_true_anomaly_deg, backward_degrees)
    % Event function to detect when spacecraft has moved backward by specified true anomaly
    kep = car2kep(x', mu);
    current_true_anomaly_deg = rad2deg(kep(6));
    
    angular_change = mod(initial_true_anomaly_deg - current_true_anomaly_deg + 180, 360) - 180;
    
    if angular_change < 0
        angular_change = angular_change + 360;
    end
    
    tolerance_deg = 0.1;
    value = backward_degrees - angular_change + tolerance_deg;
    
    isterminal = 1;
    direction = -1;
end

function positions = get_body_positions_wrapper(t_sec, t0_mjd_start, naif_ids, sp_param)
    current_mjd = t0_mjd_start + t_sec / 86400;
    num_bodies = length(naif_ids);
    positions = zeros(3, num_bodies);
    for k = 1:num_bodies
        [r, ~] = EphSS_car_spice2(naif_ids(k), current_mjd, true, sp_param);
        positions(:, k) = r(:);
    end
end