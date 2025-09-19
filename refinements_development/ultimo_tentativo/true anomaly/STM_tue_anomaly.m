function [NB_in_corrected, B_achieved_final, actual_time_diff, deltaV_magnitude, iterations_count, converged_flag, B_R_error, B_T_error, backward_duration] = STM_tue_anomaly(pars, nodein, nodeout, suppress_output)
% performSTMOptimization: Use STM to correct initial conditions for B-plane targeting
% Modified to use true anomaly-based backward propagation
%
% INPUTS:
%   pars - Parameter structure
%   nodein - Input node parameters
%   nodeout - Output node parameters  
%   pars.backward_true_anomaly_deg - Degrees to propagate backward from pericenter
%   suppress_output - (optional) If true, suppress all fprintf output (default: false)
%
% OUTPUTS:
%   NB_in_corrected - Corrected initial state vector
%   B_achieved_final - Final B-plane vector
%   actual_time_diff - Time difference at periapsis (minutes)
%   deltaV_magnitude - Required deltaV magnitude (m/s)
%   iterations_count - Number of iterations
%   converged_flag - Convergence flag (1 if converged, 0 if not)
%   B_R_error - Final B_R error (km)
%   B_T_error - Final B_T error (km)
%   backward_duration - Backward propagation duration (hours)

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
%  CALCULATE TARGET B-PLANE (from linked conics)
%  ========================================================================

[r_enceladus_at_flyby, v_enceladus_at_flyby] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, pars.INPUTS.epoch0, true, spiceParam);

[vvinfin, r0_sc_in, v0_sc_in, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodein(1), nodein(2), nodein(3), pars.INPUTS.epoch0, pars);
[vvinfout, r0_sc_out, v0_sc_out, vvga] = vinfAlphaCrank_to_VinfCARTClaudio(nodeout(1), nodeout(2), nodeout(3), pars.INPUTS.epoch0, pars);

% Calculate target B-plane parameters
e_fly = 1 + ((pars.INPUTS.Flyby.min_h + pars.Moon.EquRad) * pars.INPUTS.V_inf^2) / pars.Moon.mu;
delta_max = 2 * asin(1 / e_fly);
Rm = buildRm(r_enceladus_at_flyby, v_enceladus_at_flyby);
[vvinfouBM, delta] = infVelBeforeDefect(vvinfin, vvinfout, delta_max, vvga);

[B_vec_target, ~, B_R_target, B_T_target] = b_plane_from_asymptotes(vvinfin', vvinfouBM', pars.Moon.mu);

vvinfin_bf = (Rm * vvinfin')';
vvinfouBM_bf = (Rm * vvinfouBM')';
Energy = 0.5 * norm(vvinfin_bf)^2;
sma = - pars.Moon.mu / (2 * Energy);
ecc = 1 / (sin(delta / 2));
rp = sma * (1 - ecc);
vp = sqrt(norm(vvinfin_bf)^2 + 2 *  pars.Moon.mu / rp);
pericenter_altitude = rp - pars.Moon.EquRad;

if ~suppress_output
    fprintf('\n=== STM OPTIMIZATION FOR B-PLANE TARGETING ===\n');
    fprintf('  B_R = %.6f km   B_T = %.6f km   |B| = %.6f km\n', B_R_target, B_T_target, norm(B_vec_target));
    fprintf('  Periapse Velocity: %.3f km/s, Periapse Altitude: %.3f km\n\n', vp, pericenter_altitude);
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
NB_in = ye_bwd(1, :)';
backward_duration_sec = abs(te_bwd(1));
backward_duration_hours = backward_duration_sec / 3600;
backward_duration_days = backward_duration_sec / 86400;
initial_epoch_offset_days = pars.INPUTS.epoch0 + te_bwd(1)/86400;

% Verify final true anomaly
final_kep = car2kep(NB_in', pars.Planet.mu);
final_true_anomaly_deg = rad2deg(final_kep(6));
actual_change_deg = mod(initial_true_anomaly_deg - final_true_anomaly_deg + 180, 360) - 180;

if ~suppress_output
    fprintf('  Time duration: %.2f hours (%.3f days)\n', backward_duration_hours, backward_duration_days);
    fprintf('  Final true anomaly: %.2f degrees\n', final_true_anomaly_deg);
end

%% ========================================================================
%  ITERATIVE STM OPTIMIZATION
%  ========================================================================

% Optimization parameters
max_iterations = 25;
tolerance = 0.01; % km tolerance for B-plane
damping_factor = 0.7; % Damping for stability

% Initial guess
X_current = NB_in(1:6);
X_initial = X_current; % Store original for deltaV calculation
iteration = 0;
converged = false;

% For forward propagation, use double the backward time
flyby_tot_duration = backward_duration_sec * 2;
time_steps = pars.GroundTr.npoints;

if ~suppress_output
    fprintf('Starting STM optimization...\n');
    fprintf('Forward propagation duration: %.2f hours\n\n', flyby_tot_duration/3600);
    fprintf('Iter |  B_R Error  |  B_T Error  | Total Error |    DeltaV    | Time Diff (min)\n');
    fprintf('-----|-------------|-------------|-------------|--------------|-----------------\n');
end

while iteration < max_iterations && ~converged
    iteration = iteration + 1;
    
    % --- Propagate with STM ---
    time_vector_nbody = linspace(eps, flyby_tot_duration, time_steps)';
    
    ephem_handle = @(t_sec) get_body_positions_wrapper(t_sec, initial_epoch_offset_days, actualBodyNaifIDs, spiceParam);
    ephem_enceladus_nbody = @(t) get_body_positions_wrapper(t, initial_epoch_offset_days, pars.INPUTS.NAIFMoon, spiceParam);
    
    % No event function for this propagation - we want the full trajectory
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
    epoch_peri = initial_epoch_offset_days + time_at_peri/86400;

    time_difference_minutes = (epoch_peri - pars.INPUTS.epoch0) * 24 * 60;
    
    % --- Calculate achieved B-plane ---
    [r_enc_peri, v_enc_peri] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, epoch_peri, true, spiceParam);
    r_sc_rel = state_at_peri(1:3) - r_enc_peri';
    v_sc_rel = state_at_peri(4:6) - v_enc_peri';
    
    v_inf_sq = dot(v_sc_rel, v_sc_rel) - 2*pars.Moon.mu/norm(r_sc_rel);
    [B_R_achieved, B_T_achieved, B_vec_achieved] = b_plane_targeting(r_sc_rel, v_sc_rel, v_inf_sq, pars.Moon.mu);
    
    % --- Calculate errors ---
    delta_B_R = B_R_achieved - B_R_target;
    delta_B_T = B_T_achieved - B_T_target;
    error_norm = sqrt(delta_B_R^2 + delta_B_T^2);
    
    % Calculate cumulative deltaV
    deltaV_total = norm(X_current(4:6) - X_initial(4:6)) * 1000; % m/s
    
    if ~suppress_output
        fprintf('%3d  | %11.6f | %11.6f | %11.6f | %8.3f m/s | %14.3f\n',...
            iteration, delta_B_R, delta_B_T, error_norm, deltaV_total, time_difference_minutes);
    end
    
    % --- Check convergence ---
    if error_norm < tolerance
        converged = true;
        break;
    end
    
    % --- Compute B-plane to state sensitivity ---
    dB_dState = computeBplaneSensitivity(r_sc_rel, v_sc_rel, pars.Moon.mu);
    dB_dX0 = dB_dState * STM_at_peri;
    
    % Extract relevant sensitivities
    H = [dB_dX0(1, 4:6);   % ∂B_R/∂v₀
         dB_dX0(2, 4:6)];  % ∂B_T/∂v₀
    
    % Newton-Raphson correction
    target_error = -[delta_B_R; delta_B_T];
    delta_v = pinv(H) * target_error;
    delta_v = damping_factor * delta_v; % Apply damping
    
    % Update initial state
    X_current(4:6) = X_current(4:6) + delta_v;
end

%% ========================================================================
%  FINAL RESULTS
%  ========================================================================

if ~suppress_output
    if converged
        fprintf('\nOptimization CONVERGED in %d iterations!\n', iteration);
    else
        fprintf('\nOptimization reached maximum iterations without full convergence.\n');
    end

    % Calculate total deltaV required
    deltaV_vector = X_current(4:6) - X_initial(4:6); % km/s
    deltaV_magnitude_temp = norm(deltaV_vector) * 1000; % Convert to m/s

    fprintf('\nFinal B-plane errors:\n');
    fprintf('  B_R error: %.6f km\n', delta_B_R);
    fprintf('  B_T error: %.6f km\n', delta_B_T);
    fprintf('  Total error: %.6f km\n', error_norm);

    fprintf('\nRequired Correction:\n');
    fprintf('  DeltaV magnitude: %.3f m/s\n', deltaV_magnitude_temp);
    fprintf('  DeltaV vector: [%.6f, %.6f, %.6f] m/s\n', deltaV_vector(1)*1000, deltaV_vector(2)*1000, deltaV_vector(3)*1000);
end

% Output
NB_in_corrected = X_current;
B_achieved_final = [B_R_achieved; B_T_achieved];
actual_time_diff = time_difference_minutes; % Return the actual time difference

% Additional optional outputs for sweep analysis
deltaV_magnitude = norm(X_current(4:6) - X_initial(4:6)) * 1000; % Convert to m/s
iterations_count = iteration;
converged_flag = converged;
B_R_error = delta_B_R;
B_T_error = delta_B_T;
backward_duration = backward_duration_hours;

end

%% ========================================================================
%  TRUE ANOMALY EVENT FUNCTION
%  ========================================================================

function [value, isterminal, direction] = trueAnomalyBackwardEvent(t, x, mu, initial_true_anomaly_deg, backward_degrees)
    % Event function to detect when spacecraft has moved backward by specified true anomaly
    
    % Calculate current Keplerian elements
    kep = car2kep(x', mu);
    current_true_anomaly_deg = rad2deg(kep(6));
    
    % Calculate angular change from initial position
    angular_change = mod(initial_true_anomaly_deg - current_true_anomaly_deg + 180, 360) - 180;
    
    % For backward propagation, we expect angular_change to become positive
    % as we move backward in the orbit
    if angular_change < 0
        angular_change = angular_change + 360;
    end
    
    % Stop when we've moved backward by the desired amount (within tolerance)
    tolerance_deg = 0.1; % Stop within 0.1 degrees
    value = backward_degrees - angular_change + tolerance_deg;
    
    isterminal = 1; % Stop integration when event occurs
    direction = -1; % Detect only when approaching from positive side
end

%% ========================================================================
%  HELPER FUNCTION: Compute B-plane Sensitivity
%  ========================================================================

function dB_dState = computeBplaneSensitivity(r_rel, v_rel, mu)
    % Numerically compute ∂[B_R, B_T]/∂[r, v]
    % Returns a 2x6 matrix
    
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

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function positions = get_body_positions_wrapper(t_sec, t0_mjd_start, naif_ids, sp_param)
    current_mjd = t0_mjd_start + t_sec / 86400;
    num_bodies = length(naif_ids);
    positions = zeros(3, num_bodies);
    for k = 1:num_bodies
        [r, ~] = EphSS_car_spice2(naif_ids(k), current_mjd, true, sp_param);
        positions(:, k) = r(:);
    end
end