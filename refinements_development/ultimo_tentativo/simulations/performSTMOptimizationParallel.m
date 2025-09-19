function [NB_in_corrected, B_achieved_final, deltaV_ms, timeDiff_min, converged, iter_count] = ...
    performSTMOptimizationParallel(pars, nodein, nodeout, suppress_output)
% performSTMOptimizationParallel: Modified version for parallel execution
%
% Same as performSTMOptimization but with:
% - Optional output suppression for cleaner parallel execution
% - Additional return values for sweep analysis

if nargin < 4
    suppress_output = false;  % Default to showing output
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
vp = sqrt(norm(vvinfin_bf)^2 + 2 * pars.Moon.mu / rp);
pericenter_altitude = rp - pars.Moon.EquRad;

% Only print if not suppressed
if ~suppress_output
    fprintf('\n=== STM OPTIMIZATION FOR B-PLANE TARGETING ===\n');
    fprintf('  B_R = %.6f km   B_T = %.6f km   |B| = %.6f km\n', B_R_target, B_T_target, norm(B_vec_target));
    fprintf('  Periapse Velocity: %.3f km/s, Periapse Altitude: %.3f km\n\n', vp, pericenter_altitude);
end

%% ========================================================================
%  GET INITIAL CONDITIONS
%  ========================================================================

duration_sec = pars.t_prop_hours * 3600;
time_steps = pars.GroundTr.npoints;
time_vector_bwd = linspace(0, -duration_sec, time_steps)';

[timeLC_bwd, stateLC_bwd, ~, ~, ~] = propagateKeplerODE2(r0_sc_in', v0_sc_in', time_vector_bwd, pars.Planet.mu, []);
NB_in = stateLC_bwd(end, :)';
initial_epoch_offset_days = pars.INPUTS.epoch0 + time_vector_bwd(end)/86400;

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

% Initialize output variables
deltaV_ms = NaN;
timeDiff_min = NaN;
iter_count = 0;

if ~suppress_output
    fprintf('Starting STM optimization...\n');
    fprintf('Iter |  B_R Error  |  B_T Error  | Total Error |    DeltaV    | Time Diff (min)\n');
    fprintf('-----|-------------|-------------|-------------|--------------|-----------------\n');
end

while iteration < max_iterations && ~converged
    iteration = iteration + 1;
    
    % --- Propagate with STM ---
    flyby_tot_duration = duration_sec * 1.1;
    time_vector_nbody = linspace(eps, flyby_tot_duration, time_steps)';
    
    ephem_handle = @(t_sec) get_body_positions_wrapper(t_sec, initial_epoch_offset_days, actualBodyNaifIDs, spiceParam);
    ephem_enceladus_nbody = @(t) get_body_positions_wrapper(t, initial_epoch_offset_days, pars.INPUTS.NAIFMoon, spiceParam);
    
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
    
    % Store values for output (will be the final values when loop exits)
    deltaV_ms = deltaV_total;
    timeDiff_min = time_difference_minutes;
    iter_count = iteration;
    
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
    
    H = [dB_dX0(1, 4:6);   % ∂B_R/∂v₀
         dB_dX0(2, 4:6)];  % ∂B_T/∂v₀
    
    target_error = -[delta_B_R; delta_B_T];
    delta_v = pinv(H) * target_error;
    delta_v = damping_factor * delta_v;
    
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
    
    deltaV_vector = X_current(4:6) - X_initial(4:6);
    deltaV_magnitude = norm(deltaV_vector) * 1000;
    
    fprintf('\nFinal B-plane errors:\n');
    fprintf('  B_R error: %.6f km\n', delta_B_R);
    fprintf('  B_T error: %.6f km\n', delta_B_T);
    fprintf('  Total error: %.6f km\n', error_norm);
    
    fprintf('\nRequired Correction:\n');
    fprintf('  DeltaV magnitude: %.3f m/s\n', deltaV_magnitude);
    fprintf('  DeltaV vector: [%.6f, %.6f, %.6f] m/s\n', deltaV_vector(1)*1000, deltaV_vector(2)*1000, deltaV_vector(3)*1000);
end

% Output
NB_in_corrected = X_current;
B_achieved_final = [B_R_achieved; B_T_achieved];

end

%% ========================================================================
%  HELPER FUNCTIONS (include these as nested functions)
%  ========================================================================

function dB_dState = computeBplaneSensitivity(r_rel, v_rel, mu)
    eps_val = 1e-6;
    dB_dState = zeros(2, 6);
    
    v_inf_sq = dot(v_rel, v_rel) - 2*mu/norm(r_rel);
    [B_R_ref, B_T_ref, ~] = b_plane_targeting(r_rel, v_rel, v_inf_sq, mu);
    
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