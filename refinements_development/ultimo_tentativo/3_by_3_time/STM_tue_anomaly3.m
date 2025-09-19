function [NB_in_corrected, B_achieved_final, actual_time_diff, deltaV_magnitude, ...
          iterations_count, converged_flag, B_R_error, B_T_error, backward_duration, ...
          NB_in_state, initial_epoch_nbody, initial_epoch_lc, backward_duration_sec, ...
          mu_TBs, actualBodyNaifIDs, specialPerturbationIDs, flyby_states, ...
          T_linked_conics, propagation_duration, ephem_handles] = ...
          STM_tue_anomaly3(pars, nodein, nodeout, suppress_output)
% STM B-plane targeting with hybrid approach: B-plane first, then time correction
%
% INPUTS:
%   pars - Parameter structure
%   nodein - Input node parameters
%   nodeout - Output node parameters  
%   pars.backward_true_anomaly_deg - Degrees to propagate backward from pericenter
%   suppress_output - (optional) If true, suppress all fprintf output (default: false)
%
% OUTPUTS:
%   [Same as original function]

% Handle optional suppress_output parameter
if nargin < 4
    suppress_output = false;
end

%% ========================================================================
%  SETUP 
%  ========================================================================

% Perturbing bodies
perturbingBodyNaifIDs = pars.INPUTS.perturbingBodyNaifIDs;
specialPerturbationIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs < 0);

% Build mu_TBs list
actualBodyNaifIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs >= 0);
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

spiceParam.frame = 'J2000';
spiceParam.abcorr = 'NONE';
spiceParam.observer = '699';

%% ========================================================================
%  CALCULATE TARGET B-PLANE AND FLYBY STATES
%  ========================================================================

[r_enceladus_at_flyby, v_enceladus_at_flyby] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, pars.INPUTS.epoch0, true, spiceParam);

[vvinfin, r0_sc_in, v0_sc_in, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodein(1), nodein(2), nodein(3), pars.INPUTS.epoch0, pars);
[vvinfout, r0_sc_out, v0_sc_out, vvga] = vinfAlphaCrank_to_VinfCARTClaudio(nodeout(1), nodeout(2), nodeout(3), pars.INPUTS.epoch0, pars);

% Store flyby states for output
flyby_states.r0_sc_in = r0_sc_in;
flyby_states.v0_sc_in = v0_sc_in;
flyby_states.r0_sc_out = r0_sc_out;
flyby_states.v0_sc_out = v0_sc_out;
flyby_states.vvinfin = vvinfin;
flyby_states.vvinfout = vvinfout;

% Calculate linked conics period and propagation duration
lc_state_flyby = [r0_sc_in, v0_sc_in];
lc_kep = car2kep(lc_state_flyby, pars.Planet.mu);
T_linked_conics = 2 * pi * sqrt(lc_kep(1)^3 / pars.Planet.mu); % Period in seconds
propagation_duration = T_linked_conics/2; % Recommended propagation duration

% Calculate target B-plane parameters
e_fly = 1 + ((pars.INPUTS.Flyby.min_h + pars.Moon.EquRad) * pars.INPUTS.V_inf^2) / pars.Moon.mu;
delta_max = 2 * asin(1 / e_fly);
Rm = buildRm(r_enceladus_at_flyby, v_enceladus_at_flyby);
[vvinfouBM, delta] = infVelBeforeDefect(vvinfin, vvinfout, delta_max, vvga);

[B_vec_target, ~, B_R_target, B_T_target] = b_plane_from_asymptotes(vvinfin', vvinfouBM', pars.Moon.mu);

target_time_diff = 0; % Target is exact flyby time

vvinfin_bf = (Rm * vvinfin')';
vvinfouBM_bf = (Rm * vvinfouBM')';
Energy = 0.5 * norm(vvinfin_bf)^2;
sma = - pars.Moon.mu / (2 * Energy);
ecc = 1 / (sin(delta / 2));
rp = sma * (1 - ecc);
vp = sqrt(norm(vvinfin_bf)^2 + 2 *  pars.Moon.mu / rp);
pericenter_altitude = rp - pars.Moon.EquRad;

if ~suppress_output
    fprintf('\n=== HYBRID STM OPTIMIZATION FOR B-PLANE + TIME TARGETING ===\n');
    fprintf('  Target B_R = %.6f km   B_T = %.6f km   |B| = %.6f km\n', B_R_target, B_T_target, norm(B_vec_target));
    fprintf('  Target time error = %.1f minutes\n', target_time_diff);
    fprintf('  Periapse Velocity: %.3f km/s, Periapse Altitude: %.3f km\n', vp, pericenter_altitude);
    fprintf('  Linked Conics Period: %.2f hours\n', T_linked_conics/3600);
    fprintf('  Recommended Propagation Duration: %.2f hours\n\n', propagation_duration/3600);
end

%% ========================================================================
%  GET INITIAL CONDITIONS USING TRUE ANOMALY-BASED PROPAGATION
%  ========================================================================

% Calculate initial true anomaly at flyby conditions
initial_kep = car2kep([r0_sc_in; v0_sc_in], pars.Planet.mu);
initial_true_anomaly_deg = rad2deg(initial_kep(6));

if ~suppress_output
    fprintf('=== TRUE ANOMALY BACKWARD PROPAGATION ===\n');
    fprintf('Initial true anomaly: %.2f degrees\n', initial_true_anomaly_deg);
    fprintf('Target backward propagation: %.2f degrees\n', pars.backward_true_anomaly_deg);
end

% Create event function for true anomaly detection
event_func = @(t, x) trueAnomalyBackwardEvent(t, x, pars.Planet.mu, initial_true_anomaly_deg, pars.backward_true_anomaly_deg);

% Propagate backward with event detection
max_duration_sec = 10 * 24 * 3600; % Maximum 10 days
time_vector_bwd = linspace(0, -max_duration_sec, pars.GroundTr.npoints);

[timeLC_bwd, stateLC_bwd, te_bwd, ye_bwd, ~] = propagateKeplerODE2(r0_sc_in', v0_sc_in', time_vector_bwd, pars.Planet.mu, event_func);

if isempty(te_bwd)
    error('True anomaly event not detected within maximum propagation time');
end

% Extract results
NB_in_state = ye_bwd(1, :)'; % Original uncorrected initial state
backward_duration_sec = abs(te_bwd(1));
backward_duration_hours = backward_duration_sec / 3600;
backward_duration_days = backward_duration_sec / 86400;

% Calculate epochs
initial_epoch_lc = pars.INPUTS.epoch0; % Linked conics reference epoch (flyby)
initial_epoch_nbody = pars.INPUTS.epoch0 + te_bwd(1)/86400; % N-body initial epoch (corrected, earlier)

% Create local pars structure with computed values for potential time correction
pars_local = pars;
pars_local.initial_epoch_nbody = initial_epoch_nbody;
pars_local.mu_TBs = mu_TBs;
pars_local.actualBodyNaifIDs = actualBodyNaifIDs;
pars_local.specialPerturbationIDs = specialPerturbationIDs;

% Verify final true anomaly
final_kep = car2kep(NB_in_state', pars.Planet.mu);
final_true_anomaly_deg = rad2deg(final_kep(6));
actual_change_deg = mod(initial_true_anomaly_deg - final_true_anomaly_deg + 180, 360) - 180;

if ~suppress_output
    fprintf('  Time duration: %.2f hours (%.3f days)\n', backward_duration_hours, backward_duration_days);
    fprintf('  Final true anomaly: %.2f degrees\n', final_true_anomaly_deg);
    fprintf('  N-Body initial epoch: %.6f MJD2000\n', initial_epoch_nbody);
    fprintf('  Linked Conics reference epoch: %.6f MJD2000\n', initial_epoch_lc);
end

%% ========================================================================
%  SETUP EPHEMERIS HANDLES
%  ========================================================================

% Create ephemeris function handles for output
ephem_handles.nb_bodies = @(t_sec) get_body_positions_wrapper(t_sec, initial_epoch_nbody, actualBodyNaifIDs, spiceParam);
ephem_handles.enceladus_nb = @(t) get_body_positions_wrapper(t, initial_epoch_nbody, pars.INPUTS.NAIFMoon, spiceParam);
ephem_handles.enceladus_lc = @(t) get_body_positions_wrapper(t, initial_epoch_lc, pars.INPUTS.NAIFMoon, spiceParam);
ephem_handles.titan_lc = @(t) get_body_positions_wrapper(t, initial_epoch_lc, 606, spiceParam);

%% ========================================================================
%  HYBRID OPTIMIZATION: PHASE 1 - B-PLANE TARGETING (2x3 SYSTEM)
%  ========================================================================

% Phase 1 parameters
max_iterations_phase1 = 20;
tolerance = 0.01; % km tolerance for B-plane
damping_factor = 0.7; % Damping for stability

% Initial guess
X_current = NB_in_state(1:6);
X_initial = X_current; % Store original for deltaV calculation
iteration = 0;
phase1_converged = false;

% For forward propagation, use double the backward time
flyby_tot_duration = backward_duration_sec * 2;
time_steps = pars.GroundTr.npoints;

if ~suppress_output
    fprintf('\n=== PHASE 1: B-PLANE TARGETING (PROVEN 2x3 SYSTEM) ===\n');
    fprintf('Forward propagation duration: %.2f hours\n\n', flyby_tot_duration/3600);
    fprintf('Iter |  B_R Error  |  B_T Error  | Total Error |    DeltaV    | Time Diff (min)\n');
    fprintf('-----|-------------|-------------|-------------|--------------|-----------------\n');
end

while iteration < max_iterations_phase1 && ~phase1_converged
    iteration = iteration + 1;
    
    % --- Propagate with STM ---
    time_vector_nbody = linspace(eps, flyby_tot_duration, time_steps)';
    
    ephem_handle = ephem_handles.nb_bodies;
    ephem_enceladus_nbody = ephem_handles.enceladus_nb;
    
    [time_out, state_out, STM_history] = propagateNBodyWithSTM(X_current(1:3), X_current(4:6), ...
        time_vector_nbody, pars.Planet.mu, mu_TBs, ephem_handle, specialPerturbationIDs, []);
    
    % --- Find periapsis ---
    num_steps = length(time_out);
    dist_to_enceladus = zeros(num_steps, 1);
    for i = 1:num_steps
        r_enc_current = ephem_enceladus_nbody(time_out(i));
        dist_to_enceladus(i) = norm(state_out(i, 1:3) - r_enc_current');
    end
    
    [min_dist, idx_peri] = min(dist_to_enceladus);
    state_at_peri = state_out(idx_peri, :)';
    STM_at_peri = reshape(STM_history(idx_peri, :), 6, 6);
    time_at_peri = time_out(idx_peri);
    epoch_peri = initial_epoch_nbody + time_at_peri/86400;

    time_difference_minutes = (epoch_peri - pars.INPUTS.epoch0) * 24 * 60;
    
    % --- Calculate achieved B-plane ---
    [r_enc_peri, v_enc_peri] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, epoch_peri, true, spiceParam);
    r_sc_rel = state_at_peri(1:3) - r_enc_peri';
    v_sc_rel = state_at_peri(4:6) - v_enc_peri';
    
    v_inf_sq = dot(v_sc_rel, v_sc_rel) - 2*pars.Moon.mu/norm(r_sc_rel);
    [B_R_achieved, B_T_achieved, B_vec_achieved] = b_plane_targeting(r_sc_rel, v_sc_rel, v_inf_sq, pars.Moon.mu);
    
    % --- Calculate B-plane errors only (ignore time for Phase 1) ---
    delta_B_R = B_R_achieved - B_R_target;
    delta_B_T = B_T_achieved - B_T_target;
    error_norm = sqrt(delta_B_R^2 + delta_B_T^2);
    
    % Calculate cumulative deltaV
    deltaV_total = norm(X_current(4:6) - X_initial(4:6)) * 1000; % m/s
    
    if ~suppress_output
        fprintf('%3d  | %11.6f | %11.6f | %11.6f | %8.3f m/s | %14.3f\n',...
            iteration, delta_B_R, delta_B_T, error_norm, deltaV_total, time_difference_minutes);
    end
    
    % --- Check B-plane convergence only ---
    if abs(delta_B_R) < tolerance && abs(delta_B_T) < tolerance
        phase1_converged = true;
        break;
    end
    
    % --- Compute B-plane sensitivity (proven 2x3 system) ---
    dB_dState = computeBplaneSensitivity(r_sc_rel, v_sc_rel, pars.Moon.mu);
    dB_dX0 = dB_dState * STM_at_peri;
    
    % Extract 2x3 sensitivity matrix (B_R, B_T vs v_x, v_y, v_z)
    H = dB_dX0(:, 4:6);  % 2x3 matrix
    
    % 2x3 system: use pseudo-inverse (proven approach)
    target_error = -[delta_B_R; delta_B_T];
    delta_v = pinv(H) * target_error;
    delta_v = damping_factor * delta_v; % Apply damping
    
    % Update initial state
    X_current(4:6) = X_current(4:6) + delta_v;
end

% Store Phase 1 results
X_after_phase1 = X_current;
phase1_time_error = time_difference_minutes - target_time_diff;
phase1_deltaV = norm(X_current(4:6) - X_initial(4:6)) * 1000;
phase1_B_R_error = delta_B_R;
phase1_B_T_error = delta_B_T;

if ~suppress_output
    fprintf('\n=== PHASE 1 RESULTS ===\n');
    if phase1_converged
        fprintf('B-plane targeting CONVERGED in %d iterations\n', iteration);
    else
        fprintf('B-plane targeting reached maximum iterations\n');
    end
    fprintf('Final B-plane errors: B_R = %.6f km, B_T = %.6f km\n', delta_B_R, delta_B_T);
    fprintf('Time error after B-plane correction: %.3f minutes\n', phase1_time_error);
    fprintf('Phase 1 DeltaV: %.3f m/s\n', phase1_deltaV);
end

%% ========================================================================
%  HYBRID OPTIMIZATION: PHASE 2 - TIME CORRECTION (IF NEEDED)
%  ========================================================================

% Phase 2 parameters
time_tolerance = 5.0; % minutes - acceptable time error
max_time_correction = 0.005; % km/s - maximum velocity correction for time
phase2_applied = false;
phase2_converged = false;
final_time_error = phase1_time_error;

% Only apply Phase 2 if:
% 1. Phase 1 converged (B-plane targeting successful)
% 2. Time error is larger than tolerance
% 3. Time error is reasonable (not huge - indicates other problems)
if phase1_converged && abs(phase1_time_error) > time_tolerance && abs(phase1_time_error) < 300
    if ~suppress_output
        fprintf('\n=== PHASE 2: TIME CORRECTION ===\n');
        fprintf('Attempting to reduce time error from %.3f to ≤ %.1f minutes\n', phase1_time_error, time_tolerance);
    end
    
    % Simple time sensitivity analysis
    % Test each velocity component to find which affects timing most
    eps_test = 1e-6; % Small test perturbation (km/s)
    time_sensitivities = zeros(1, 3);
    
    for comp = 1:3
        X_test = X_current;
        X_test(3 + comp) = X_test(3 + comp) + eps_test;
        
        % Quick propagation to check time sensitivity
        try
            [~, ~, time_test] = computePerturbedBplaneAndTime(X_test, pars_local, ephem_handles);
            time_sensitivities(comp) = (time_test - time_difference_minutes) / eps_test; % minutes per km/s
        catch
            time_sensitivities(comp) = 0;
        end
    end
    
    % Find component with strongest time sensitivity
    [max_sensitivity, best_component] = max(abs(time_sensitivities));
    
    if max_sensitivity > 1e-6 % Sensitivity is meaningful
        % Calculate required velocity correction
        required_time_change = -phase1_time_error; % minutes (negative to correct the error)
        velocity_correction = required_time_change / time_sensitivities(best_component); % km/s
        
        % Limit correction to reasonable magnitude
        velocity_correction = sign(velocity_correction) * min(abs(velocity_correction), max_time_correction);
        
        if ~suppress_output
            fprintf('Strongest time sensitivity: %.2f min/(km/s) on v_%d\n', time_sensitivities(best_component), best_component);
            fprintf('Applying %.6f km/s correction to v_%d\n', velocity_correction, best_component);
        end
        
        % Apply correction
        X_current(3 + best_component) = X_current(3 + best_component) + velocity_correction;
        
        % Verify the correction by re-propagating
        [time_out_verify, state_out_verify, ~] = propagateNBodyWithSTM(X_current(1:3), X_current(4:6), ...
            time_vector_nbody, pars.Planet.mu, mu_TBs, ephem_handle, specialPerturbationIDs, []);
        
        % Find new periapsis timing
        dist_to_enceladus_verify = zeros(length(time_out_verify), 1);
        for i = 1:length(time_out_verify)
            r_enc_current = ephem_enceladus_nbody(time_out_verify(i));
            dist_to_enceladus_verify(i) = norm(state_out_verify(i, 1:3) - r_enc_current');
        end
        
        [~, idx_peri_verify] = min(dist_to_enceladus_verify);
        state_at_peri_verify = state_out_verify(idx_peri_verify, :)';
        time_at_peri_verify = time_out_verify(idx_peri_verify);
        epoch_peri_verify = initial_epoch_nbody + time_at_peri_verify/86400;
        new_time_difference = (epoch_peri_verify - pars.INPUTS.epoch0) * 24 * 60;
        final_time_error = new_time_difference - target_time_diff;
        
        % Check B-plane degradation after time correction
        [r_enc_peri_verify, v_enc_peri_verify] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, epoch_peri_verify, true, spiceParam);
        r_sc_rel_verify = state_at_peri_verify(1:3) - r_enc_peri_verify';
        v_sc_rel_verify = state_at_peri_verify(4:6) - v_enc_peri_verify';
        v_inf_sq_verify = dot(v_sc_rel_verify, v_sc_rel_verify) - 2*pars.Moon.mu/norm(r_sc_rel_verify);
        [B_R_verify, B_T_verify, ~] = b_plane_targeting(r_sc_rel_verify, v_sc_rel_verify, v_inf_sq_verify, pars.Moon.mu);
        
        B_R_error_after = B_R_verify - B_R_target;
        B_T_error_after = B_T_verify - B_T_target;
        bplane_degradation = sqrt(B_R_error_after^2 + B_T_error_after^2) - sqrt(phase1_B_R_error^2 + phase1_B_T_error^2);
        
        % Accept correction if it improves time without significantly degrading B-plane
        time_improvement = abs(final_time_error) < abs(phase1_time_error);
        bplane_acceptable = abs(B_R_error_after) < tolerance && abs(B_T_error_after) < tolerance;
        
        if time_improvement && bplane_acceptable
            phase2_applied = true;
            time_difference_minutes = new_time_difference;
            delta_B_R = B_R_error_after;  % Update for final output
            delta_B_T = B_T_error_after;
            B_R_achieved = B_R_verify;   % Update for final output  
            B_T_achieved = B_T_verify;
            
            if abs(final_time_error) <= time_tolerance
                phase2_converged = true;
            end
            
            if ~suppress_output
                fprintf('Time correction SUCCESSFUL:\n');
                fprintf('  Time error: %.3f → %.3f minutes\n', phase1_time_error, final_time_error);
                fprintf('  B-plane degradation: %.6f km\n', bplane_degradation);
                fprintf('  Final B-plane errors: B_R = %.6f km, B_T = %.6f km\n', B_R_error_after, B_T_error_after);
            end
        else
            % Revert if correction made things worse
            X_current = X_after_phase1;
            final_time_error = phase1_time_error;
            if ~suppress_output
                fprintf('Time correction REJECTED (time_improved=%d, bplane_ok=%d)\n', time_improvement, bplane_acceptable);
                fprintf('  Reverted to Phase 1 solution\n');
            end
        end
    else
        if ~suppress_output
            fprintf('Time sensitivity too weak (%.2e min/(km/s)), skipping correction\n', max_sensitivity);
        end
    end
else
    if ~suppress_output
        if ~phase1_converged
            fprintf('Phase 2 skipped: Phase 1 did not converge\n');
        elseif abs(phase1_time_error) <= time_tolerance
            fprintf('Phase 2 not needed: time error (%.3f min) within tolerance\n', phase1_time_error);
        else
            fprintf('Phase 2 skipped: time error too large (%.3f min), indicates other issues\n', phase1_time_error);
        end
    end
end

%% ========================================================================
%  FINAL RESULTS AND OUTPUT
%  ========================================================================

% Determine overall convergence
converged = phase1_converged && (phase2_converged || abs(final_time_error) <= time_tolerance);
total_deltaV = norm(X_current(4:6) - X_initial(4:6)) * 1000;
total_iterations = iteration + phase2_applied;

if ~suppress_output
    fprintf('\n=== HYBRID OPTIMIZATION SUMMARY ===\n');
    fprintf('Phase 1 (B-plane): %s in %d iterations\n', ...
            ternary(phase1_converged, 'CONVERGED', 'Max iterations'), iteration);
    fprintf('Phase 2 (time): %s\n', ...
            ternary(~phase1_converged || abs(phase1_time_error) <= time_tolerance, 'NOT NEEDED', ...
                   ternary(phase2_applied, ...
                          ternary(phase2_converged, 'CONVERGED', 'IMPROVED'), 'FAILED')));
    fprintf('Overall status: %s\n', ternary(converged, 'CONVERGED', 'PARTIAL'));
    fprintf('Total iterations: %d\n', total_iterations);
    fprintf('Final errors:\n');
    fprintf('  B_R = %.6f km, B_T = %.6f km\n', delta_B_R, delta_B_T);
    fprintf('  Time = %.3f minutes\n', final_time_error);
    fprintf('Total DeltaV: %.3f m/s\n', total_deltaV);
    
    % Performance comparison note
    fprintf('\nNote: This hybrid approach prioritizes B-plane precision with optional time improvement.\n');
    fprintf('Pure 2x3 system typically achieves B-plane errors < 0.01 km with 70-100 m/s DeltaV.\n');
end

% Output assignments
NB_in_corrected = X_current;
B_achieved_final = [B_R_achieved; B_T_achieved];
actual_time_diff = final_time_error; % Return final time error (could be phase1 or phase2)
deltaV_magnitude = total_deltaV;
iterations_count = total_iterations;
converged_flag = converged;
B_R_error = delta_B_R;
B_T_error = delta_B_T;
backward_duration = backward_duration_hours;

end

%% ========================================================================
%  SUPPORT FUNCTIONS
%  ========================================================================

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

function dB_dState = computeBplaneSensitivity(r_rel, v_rel, mu)
    % Numerically compute ∂[B_R, B_T]/∂[r, v] - proven method
    eps_val = 1e-6; % Perturbation size
    dB_dState = zeros(2, 6);
    
    % Reference B-plane
    v_inf_sq = dot(v_rel, v_rel) - 2*mu/norm(r_rel);
    [B_R_ref, B_T_ref, ~] = b_plane_targeting(r_rel, v_rel, v_inf_sq, mu);
    
    % Perturb each state component
    for i = 1:6
        state_pert = [r_rel; v_rel];
        state_pert(i) = state_pert(i) + eps_val;
        
        r_pert = state_pert(1:3);
        v_pert = state_pert(4:6);
        
        v_inf_sq_pert = dot(v_pert, v_pert) - 2*mu/norm(r_pert);
        [B_R_pert, B_T_pert, ~] = b_plane_targeting(r_pert, v_pert, v_inf_sq_pert, mu);
        
        dB_dState(1, i) = (B_R_pert - B_R_ref) / eps_val;
        dB_dState(2, i) = (B_T_pert - B_T_ref) / eps_val;
    end
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

function [B_R_achieved, B_T_achieved, time_diff_minutes] = computePerturbedBplaneAndTime(X_perturbed, pars, ephem_handles)
    % Simplified propagation for time sensitivity computation
    flyby_duration = 24 * 3600; % 1 day propagation for speed
    time_steps = 300; % Reduced points for speed
    time_vector = linspace(eps, flyby_duration, time_steps)';
    
    [time_out, state_out, ~] = propagateNBodyWithSTM(X_perturbed(1:3), X_perturbed(4:6), ...
        time_vector, pars.Planet.mu, pars.mu_TBs, ephem_handles.nb_bodies, pars.specialPerturbationIDs, []);
    
    % Find pericenter
    dist_to_enceladus = zeros(length(time_out), 1);
    for i = 1:length(time_out)
        r_enc_current = ephem_handles.enceladus_nb(time_out(i));
        dist_to_enceladus(i) = norm(state_out(i, 1:3) - r_enc_current');
    end
    
    [~, idx_peri] = min(dist_to_enceladus);
    state_at_peri = state_out(idx_peri, :)';
    time_at_peri = time_out(idx_peri);
    
    epoch_peri = pars.initial_epoch_nbody + time_at_peri/86400;
    time_diff_minutes = (epoch_peri - pars.INPUTS.epoch0) * 24 * 60;
    
    % Calculate B-plane parameters
    spiceParam.frame = 'J2000';
    spiceParam.abcorr = 'NONE';
    spiceParam.observer = '699';
    
    [r_enc_peri, v_enc_peri] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, epoch_peri, true, spiceParam);
    r_sc_rel = state_at_peri(1:3) - r_enc_peri';
    v_sc_rel = state_at_peri(4:6) - v_enc_peri';
    
    v_inf_sq = dot(v_sc_rel, v_sc_rel) - 2*pars.Moon.mu/norm(r_sc_rel);
    [B_R_achieved, B_T_achieved, ~] = b_plane_targeting(r_sc_rel, v_sc_rel, v_inf_sq, pars.Moon.mu);
end

function result = ternary(condition, true_val, false_val)
    % Simple ternary operator implementation
    if condition
        result = true_val;
    else
        result = false_val;
    end
end