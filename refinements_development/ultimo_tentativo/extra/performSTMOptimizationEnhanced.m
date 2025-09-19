function [NB_in_corrected, B_achieved_final] = performSTMOptimizationEnhanced(pars, nodein, nodeout)
    % Enhanced version with better handling of close encounters
    %
    % This version includes:
    % - DeltaV reporting
    % - Smart handling of Enceladus when it's a perturber
    % - Better convergence diagnostics
    
    %% ========================================================================
    %  SETUP
    %  ========================================================================
    
    soi_multiplier = 64;
    pars.t_prop_hours = 240;
    pars.GroundTr.npoints = 30e3;
    
    % Identify if Enceladus is among the perturbers
    perturbingBodyNaifIDs = pars.INPUTS.perturbingBodyNaifIDs;
    has_enceladus_perturber = any(perturbingBodyNaifIDs == 602);
    
    if has_enceladus_perturber
        fprintf('\nWARNING: Enceladus (602) detected as a perturber.\n');
        fprintf('         Close encounter singularities will be handled automatically.\n\n');
    end
    
    % Separate special and actual body perturbations
    specialPerturbationIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs < 0);
    actualBodyNaifIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs >= 0);
    
    % Build mu_TBs list with body identification
    mu_TBs = zeros(1, length(actualBodyNaifIDs));
    body_names = cell(1, length(actualBodyNaifIDs));
    
    for i = 1:length(actualBodyNaifIDs)
        id = actualBodyNaifIDs(i);
        switch id
            case 10
                mu_TBs(i) = getAstroConstants('Sun', 'Mu');
                body_names{i} = 'Sun';
            case 5
                mu_TBs(i) = getAstroConstants('Jupiter', 'Mu');
                body_names{i} = 'Jupiter';
            case 601
                [~, mu, ~, ~] = satMoonsConstants(0); mu_TBs(i) = mu;
                body_names{i} = 'Mimas';
            case 602
                [~, mu, ~, ~] = satMoonsConstants(1); mu_TBs(i) = mu;
                body_names{i} = 'Enceladus';
            case 603
                [~, mu, ~, ~] = satMoonsConstants(2); mu_TBs(i) = mu;
                body_names{i} = 'Tethys';
            case 604
                [~, mu, ~, ~] = satMoonsConstants(3); mu_TBs(i) = mu;
                body_names{i} = 'Dione';
            case 605
                [~, mu, ~, ~] = satMoonsConstants(4); mu_TBs(i) = mu;
                body_names{i} = 'Rhea';
            case 606
                [~, mu, ~, ~] = satMoonsConstants(5); mu_TBs(i) = mu;
                body_names{i} = 'Titan';
            case 607
                mu_TBs(i) = 0.374;
                body_names{i} = 'Hyperion';
            case 608
                mu_TBs(i) = 120.4;
                body_names{i} = 'Iapetus';
            otherwise
                mu_TBs(i) = 0;
                body_names{i} = sprintf('Body_%d', id);
        end
    end
    
    % Print perturber configuration
    fprintf('Perturber Configuration:\n');
    fprintf('------------------------\n');
    if ~isempty(specialPerturbationIDs)
        fprintf('Zonal harmonics: ');
        for j = specialPerturbationIDs
            fprintf('J%d ', abs(j));
        end
        fprintf('\n');
    end
    if ~isempty(actualBodyNaifIDs)
        fprintf('N-body perturbers:\n');
        for i = 1:length(actualBodyNaifIDs)
            fprintf('  - %s (ID: %d, mu: %.3f km³/s²)\n', body_names{i}, actualBodyNaifIDs(i), mu_TBs(i));
        end
    end
    fprintf('\n');
    
    spiceParam.frame = 'J2000';
    spiceParam.abcorr = 'NONE';
    spiceParam.observer = '699';
    
    %% ========================================================================
    %  CALCULATE TARGET B-PLANE
    %  ========================================================================
    
    [r_enceladus_at_flyby, v_enceladus_at_flyby] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, pars.INPUTS.epoch0, true, spiceParam);
    
    [vvinfin, r0_sc_in, v0_sc_in, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodein(1), nodein(2), nodein(3), pars.INPUTS.epoch0, pars);
    [vvinfout, r0_sc_out, v0_sc_out, vvga] = vinfAlphaCrank_to_VinfCARTClaudio(nodeout(1), nodeout(2), nodeout(3), pars.INPUTS.epoch0, pars);
    
    % Calculate target B-plane
    e_fly = 1 + ((pars.INPUTS.Flyby.min_h + pars.Moon.EquRad) * pars.INPUTS.V_inf^2) / pars.Moon.mu;
    delta_max = 2 * asin(1 / e_fly);
    Rm = buildRm(r_enceladus_at_flyby, v_enceladus_at_flyby);
    [vvinfouBM, delta] = infVelBeforeDefect(vvinfin, vvinfout, delta_max, vvga);
    vvinfin_bf = (Rm * vvinfin')';
    vvinfouBM_bf = (Rm * vvinfouBM')';
    
    [B_vec_target, ~, B_R_target, B_T_target] = b_plane_from_asymptotes(vvinfin', vvinfouBM', pars.Moon.mu);
    
    fprintf('=== STM OPTIMIZATION FOR B-PLANE TARGETING ===\n');
    fprintf('Target B-plane parameters:\n');
    fprintf('  B_R = %.6f km\n', B_R_target);
    fprintf('  B_T = %.6f km\n', B_T_target);
    fprintf('  |B| = %.6f km\n\n', norm(B_vec_target));
    
    %% ========================================================================
    %  GET INITIAL CONDITIONS AT -64 SOI
    %  ========================================================================
    
    duration_sec = pars.t_prop_hours * 3600;
    time_steps = pars.GroundTr.npoints;
    time_vector_bwd = linspace(0, -duration_sec, time_steps)';
    
    soi_64_dist = pars.INPUTS.maxPropagationDistance;
    ephem_enceladus = @(t) get_body_positions_wrapper(t, pars.INPUTS.epoch0, pars.INPUTS.NAIFMoon, spiceParam);
    soi_event_func = @(t, x) soiCrossingEvent(t, x, ephem_enceladus, soi_64_dist);
    
    [~, ~, te_bwd, ye_bwd, ~] = propagateKeplerODE2(r0_sc_in', v0_sc_in', time_vector_bwd, pars.Planet.mu, soi_event_func);
    
    NB_in = ye_bwd(1, :)';
    initial_epoch_offset_days = pars.INPUTS.epoch0 + te_bwd/86400;
    
    %% ========================================================================
    %  ITERATIVE STM OPTIMIZATION
    %  ========================================================================
    
    % Optimization parameters
    max_iterations = 20;
    tolerance = 1; % km
    damping_factor = 0.7;
    min_damping = 0.3;
    max_damping = 1.0;
    
    % Initial state
    X_current = NB_in(1:6);
    X_initial = X_current;
    iteration = 0;
    converged = false;
    error_history = [];
    
    fprintf('Starting STM optimization...\n');
    fprintf('Iter |  B_R Error  |  B_T Error  |  Total Error |   DeltaV    | Damping\n');
    fprintf('-----|-------------|-------------|--------------|-------------|--------\n');
    
    while iteration < max_iterations && ~converged
        iteration = iteration + 1;
        
        % Propagate with STM
        flyby_duration = duration_sec * 2;
        time_vector_nbody = linspace(eps, flyby_duration, time_steps)';
        
        % Create enhanced ephemeris handle that tracks close encounters
        ephem_handle = @(t_sec) get_body_positions_wrapper(t_sec, initial_epoch_offset_days, actualBodyNaifIDs, spiceParam);
        ephem_enceladus_nbody = @(t) get_body_positions_wrapper(t, initial_epoch_offset_days, pars.INPUTS.NAIFMoon, spiceParam);
        
        % Propagate with STM (using the smart singularity avoidance)
        try
            [time_out, state_out, STM_history] = propagateNBodyWithSTMSmart(...
                X_current(1:3), X_current(4:6), time_vector_nbody, ...
                pars.Planet.mu, mu_TBs, ephem_handle, specialPerturbationIDs, ...
                actualBodyNaifIDs, body_names);
        catch ME
            fprintf('\nERROR during propagation at iteration %d:\n', iteration);
            fprintf('  %s\n', ME.message);
            fprintf('  Trying with reduced time steps...\n');
            
            % Try with more time steps for better numerical stability
            time_vector_nbody = linspace(eps, flyby_duration, time_steps*2)';
            [time_out, state_out, STM_history] = propagateNBodyWithSTMSmart(...
                X_current(1:3), X_current(4:6), time_vector_nbody, ...
                pars.Planet.mu, mu_TBs, ephem_handle, specialPerturbationIDs, ...
                actualBodyNaifIDs, body_names);
        end
        
        % Find periapsis
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
        
        % Calculate achieved B-plane
        [r_enc_peri, v_enc_peri] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, epoch_peri, true, spiceParam);
        r_sc_rel = state_at_peri(1:3) - r_enc_peri';
        v_sc_rel = state_at_peri(4:6) - v_enc_peri';
        
        v_inf_sq = dot(v_sc_rel, v_sc_rel) - 2*pars.Moon.mu/norm(r_sc_rel);
        
        % Check for valid v_inf_sq
        if v_inf_sq <= 0
            fprintf('\nWARNING: Negative v_inf^2 detected. Spacecraft may be captured.\n');
            fprintf('         Adjusting to use minimum positive value.\n');
            v_inf_sq = 1e-6;
        end
        
        [B_R_achieved, B_T_achieved, B_vec_achieved] = b_plane_targeting(r_sc_rel, v_sc_rel, v_inf_sq, pars.Moon.mu);
        
        % Calculate errors
        delta_B_R = B_R_achieved - B_R_target;
        delta_B_T = B_T_achieved - B_T_target;
        error_norm = sqrt(delta_B_R^2 + delta_B_T^2);
        error_history(iteration) = error_norm;
        
        % Calculate deltaV
        deltaV_total = norm(X_current(4:6) - X_initial(4:6)) * 1000; % m/s
        
        % Adaptive damping based on convergence behavior
        if iteration > 1
            if error_norm > error_history(iteration-1)
                % Error increased, reduce damping
                damping_factor = max(min_damping, damping_factor * 0.7);
            elseif error_norm < 0.5 * error_history(iteration-1)
                % Good convergence, increase damping
                damping_factor = min(max_damping, damping_factor * 1.2);
            end
        end
        
        fprintf('%3d  | %11.6f | %11.6f | %11.6f | %8.3f m/s | %.3f\n', ...
            iteration, delta_B_R, delta_B_T, error_norm, deltaV_total, damping_factor);
        
        % Check convergence
        if error_norm < tolerance
            converged = true;
            break;
        end
        
        % Check for stagnation
        if iteration > 3
            recent_errors = error_history(max(1,iteration-2):iteration);
            if std(recent_errors) < tolerance/10
                fprintf('\nOptimization appears to be stagnating. Stopping.\n');
                break;
            end
        end
        
        % Compute B-plane sensitivity
        dB_dState = computeBplaneSensitivity(r_sc_rel, v_sc_rel, pars.Moon.mu);
        
        % Check STM condition number
        cond_STM = cond(STM_at_peri);
        if cond_STM > 1e10
            fprintf('WARNING: STM condition number is very large (%.2e). Results may be inaccurate.\n', cond_STM);
        end
        
        % Map to initial conditions
        dB_dX0 = dB_dState * STM_at_peri;
        
        % Extract velocity sensitivities
        H = [dB_dX0(1, 4:6); dB_dX0(2, 4:6)];
        
        % Check sensitivity matrix condition
        cond_H = cond(H);
        if cond_H > 1e8
            fprintf('WARNING: Sensitivity matrix is ill-conditioned (cond = %.2e).\n', cond_H);
            fprintf('         Using regularization...\n');
            % Add regularization
            H = H + 1e-8 * eye(size(H));
        end
        
        % Newton-Raphson correction
        target_error = -[delta_B_R; delta_B_T];
        delta_v = pinv(H) * target_error;
        
        % Limit maximum correction
        max_dv = 0.1; % km/s maximum per iteration
        if norm(delta_v) > max_dv
            delta_v = delta_v * max_dv / norm(delta_v);
            fprintf('  Note: Correction limited to %.1f m/s\n', max_dv * 1000);
        end
        
        delta_v = damping_factor * delta_v;
        
        % Update state
        X_current(4:6) = X_current(4:6) + delta_v;
    end
    
    %% ========================================================================
    %  FINAL RESULTS
    %  ========================================================================
    
    if converged
        fprintf('\n=== OPTIMIZATION CONVERGED in %d iterations! ===\n', iteration);
    else
        fprintf('\n=== Optimization stopped after %d iterations ===\n', iteration);
    end
    
    % Calculate final deltaV
    deltaV_vector = X_current(4:6) - X_initial(4:6);
    deltaV_magnitude = norm(deltaV_vector) * 1000; % m/s
    
    fprintf('\nFinal B-plane Achievement:\n');
    fprintf('  B_R error: %.6f km (%.2f%%)\n', delta_B_R, abs(delta_B_R/B_R_target)*100);
    fprintf('  B_T error: %.6f km (%.2f%%)\n', delta_B_T, abs(delta_B_T/B_T_target)*100);
    fprintf('  Total error: %.6f km\n', error_norm);
    
    fprintf('\nRequired Trajectory Correction:\n');
    fprintf('  DeltaV magnitude: %.3f m/s\n', deltaV_magnitude);
    fprintf('  DeltaV components: [%.3f, %.3f, %.3f] m/s\n', ...
        deltaV_vector(1)*1000, deltaV_vector(2)*1000, deltaV_vector(3)*1000);
    fprintf('  Applied at: %.1f Enceladus SOI before encounter\n', soi_multiplier);
    fprintf('  Periapsis altitude achieved: %.3f km\n', min_dist - pars.Moon.EquRad);
    
    % Outputs
    NB_in_corrected = X_current;
    B_achieved_final = [B_R_achieved; B_T_achieved];

    %% ========================================================================
    %  NESTED HELPER FUNCTIONS START HERE
    %  ========================================================================
    
    function [tt, yy, STM_history] = propagateNBodyWithSTMSmart(rvec, vvec, timevector, ...
        muCB, mu_TBs, ephem_func, specialPerturbationIDs, bodyIDs, bodyNames)
        
        % Initial STM
        STM0 = eye(6);
        STM0_vec = reshape(STM0, 36, 1);
        x0_augmented = [rvec; vvec; STM0_vec];
        
        % Enhanced EOM with diagnostics
        close_encounter_logged = false;
        F = @(t, x) nBodyEOMWithSTMSmart(t, x, muCB, mu_TBs, ephem_func, ...
            specialPerturbationIDs, bodyIDs, bodyNames, close_encounter_logged);
        
        % ODE options with better tolerances for close encounters
        options = odeset('RelTol', 1e-12, 'AbsTol', 1e-14, 'MaxStep', 100);
        
        % Propagate
        [tt, xx_augmented] = ode45(F, timevector, x0_augmented, options);
        
        % Extract results
        yy = xx_augmented(:, 1:6);
        STM_history = xx_augmented(:, 7:42);
    end
    
    function dxdt = nBodyEOMWithSTMSmart(t, x, muCB, mu_TBs, ephem_func, ...
        specialPerturbationIDs, bodyIDs, bodyNames, close_encounter_logged)
        
        r_sc = x(1:3);
        v_sc = x(4:6);
        STM_vec = x(7:42);
        STM = reshape(STM_vec, 6, 6);
        
        % Central body acceleration
        a_CB = -muCB * r_sc / norm(r_sc)^3;
        
        % Zonal harmonics
        a_zonal_pert = [0; 0; 0];
        if ~isempty(specialPerturbationIDs)
            j_term_types = abs(specialPerturbationIDs);
            % NOTE: You will need to provide your own 'calculate_zonal_perturbations' function
            % if it's not already on your MATLAB path.
            all_zonal_accels = calculate_zonal_perturbations(r_sc, j_term_types);
            a_zonal_pert = sum(all_zonal_accels, 2);
        end
        
        % N-body with smart handling
        a_nbody_pert = [0; 0; 0];
        min_dist_thresh = 50; % km - increased threshold for Enceladus
        
        if ~isempty(mu_TBs)
            r_TBs = ephem_func(t);
            for i = 1:length(mu_TBs)
                mu_tb = mu_TBs(i);
                r_tb = r_TBs(:, i);
                rtb2sc = r_sc - r_tb;
                dist = norm(rtb2sc);
                
                if bodyIDs(i) == 602
                    if dist < 255 % enecleadus equatorial radius
                        continue; % Skip Enceladus perturbation when very close
                    end
                elseif dist < min_dist_thresh
                    continue; % Skip other bodies when too close
                end
                
                % Normal perturbation calculation
                pert_i = -mu_tb * (rtb2sc / dist^3) + mu_tb * (r_tb / norm(r_tb)^3);
                a_nbody_pert = a_nbody_pert + pert_i;
            end
        end
        
        a_total = a_CB + a_zonal_pert + a_nbody_pert;
        
        % Compute A matrix with same smart handling
        A = computeAMatrixSmart(r_sc, muCB, mu_TBs, ephem_func, t, ...
            specialPerturbationIDs, bodyIDs);
        
        % STM evolution
        STM_dot = A * STM;
        STM_dot_vec = reshape(STM_dot, 36, 1);
        
        dxdt = [v_sc; a_total; STM_dot_vec];
    end
    
    function A = computeAMatrixSmart(r, muCB, mu_TBs, ephem_func, t, ...
        specialPerturbationIDs, bodyIDs)
        
        A = zeros(6, 6);
        A(1:3, 4:6) = eye(3);
        
        % Central body gradient
        r_norm = norm(r);
        I3 = eye(3);
        rr_dyad = r * r';
        G_CB = -muCB/r_norm^3 * I3 + 3*muCB/r_norm^5 * rr_dyad;
        
        % J2 gradient
        G_J2 = zeros(3, 3);
        if ~isempty(specialPerturbationIDs) && any(abs(specialPerturbationIDs) == 2)
            G_J2 = computeJ2GravityGradient(r, muCB);
        end
        
        % N-body gradients with smart handling
        G_nbody = zeros(3, 3);
        min_dist_thresh = 50;
        
        if ~isempty(mu_TBs)
            r_TBs = ephem_func(t);
            for i = 1:length(mu_TBs)
                mu_tb = mu_TBs(i);
                r_tb = r_TBs(:, i);
                r_rel = r - r_tb;
                dist = norm(r_rel);
                
                if bodyIDs(i) == 602 && dist < min_dist_thresh * 2
                    continue;
                elseif dist < min_dist_thresh
                    continue;
                end
                
                G_tb = -mu_tb/dist^3 * I3 + 3*mu_tb/dist^5 * (r_rel * r_rel');
                G_nbody = G_nbody + G_tb;
            end
        end
        
        A(4:6, 1:3) = G_CB + G_J2 + G_nbody;
    end

    % ****** MOVED FUNCTION IS HERE ******
    function G_J2 = computeJ2GravityGradient(r, muCB)
        % Compute J2 gravity gradient contribution
        % This is simplified - you should use your actual J2 model
        
        % Saturn J2 parameters
        J2 = 16298e-6;  % Saturn J2
        R_ref = 60268;  % Saturn reference radius (km)
        
        x = r(1); y = r(2); z = r(3);
        r_norm = norm(r);
        r2 = r_norm^2;
        r5 = r_norm^5;
        r7 = r_norm^7;
        
        % J2 gravity gradient terms (simplified)
        factor = 3 * muCB * J2 * R_ref^2 / 2;
        
        G_J2 = zeros(3, 3);
        
        % Diagonal terms
        G_J2(1,1) = factor * (3/r5 - 30*x^2/r7 + 35*z^2*x^2/r7^2);
        G_J2(2,2) = factor * (3/r5 - 30*y^2/r7 + 35*z^2*y^2/r7^2);
        G_J2(3,3) = factor * (-9/r5 + 30*z^2/r7 - 35*z^4/r7^2);
        
        % Off-diagonal terms
        G_J2(1,2) = factor * (-30*x*y/r7 + 35*z^2*x*y/r7^2);
        G_J2(2,1) = G_J2(1,2);
        G_J2(1,3) = factor * (30*x*z/r7 - 35*z^3*x/r7^2);
        G_J2(3,1) = G_J2(1,3);
        G_J2(2,3) = factor * (30*y*z/r7 - 35*z^3*y/r7^2);
        G_J2(3,2) = G_J2(2,3);
    end

function [value, isterminal, direction] = soiCrossingEvent(t, x, ephem_enceladus_handle, target_distance)
    r_sc = x(1:3);
    r_enceladus = ephem_enceladus_handle(t);
    current_distance = norm(r_sc - r_enceladus);
    value = current_distance - target_distance;
    isterminal = 1;
    direction = 0;
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

end % This is the end of the main performSTMOptimizationEnhanced function