%% STM B-Plane Targeting with Three-Trajectory Animation
clc;
clear all;
close all;

%% ========================================================================
% 1. SETUP
% ========================================================================
pars.backward_true_anomaly_deg = 150; % Degrees to propagate backward from pericenter
pars.GroundTr.npoints = 1000; % More points for smoother animation
% pars.INPUTS.perturbingBodyNaifIDs = [-2,10];
% pars.INPUTS.perturbingBodyNaifIDs = [-2,10,602, 606]; % J2, Enceladus, Sun, Titan
pars.INPUTS.perturbingBodyNaifIDs = [-2,10, 606]; % J2, Sun, Enceladus, Titan

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
% nodein = [4, deg2rad(8.6918), deg2rad(-86.9406)];
% nodeout = [4, deg2rad(8.6918), deg2rad(-88.1610)];

% Flyby nodes - Partial-COT 1 (O/I) - 1st Flyby
nodein = [4, 0.15, deg2rad(0)];
nodeout = [4, 0.15, deg2rad(1)];

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
