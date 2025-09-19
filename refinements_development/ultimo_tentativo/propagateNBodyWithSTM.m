function [tt, yy, STM_history, te, ye, ie] = propagateNBodyWithSTM(rvec, vvec, timevector, muCB, mu_TBs, ephem_func, specialPerturbationIDs, event_func_handle)
% propagateNBodyWithSTM: N-body propagator with State Transition Matrix
%
% This function propagates spacecraft dynamics including STM for sensitivity analysis
%
% INPUTS:
%   rvec        - [3x1] Initial position vector (km)
%   vvec        - [3x1] Initial velocity vector (km/s)
%   timevector  - [Nx1] Time vector (s)
%   muCB        - Gravitational parameter of central body (km^3/s^2)
%   mu_TBs      - [Mx1] Gravitational parameters of perturbing bodies
%   ephem_func  - Function handle for ephemerides
%   specialPerturbationIDs - IDs for special perturbations (J2, etc.)
%   event_func_handle - Optional event detection function
%
% OUTPUTS:
%   tt          - Output time vector
%   yy          - [Nx6] State history [r; v]
%   STM_history - [Nx36] STM history (reshaped as column vector)
%   te, ye, ie  - Event detection outputs

if nargin < 7
    specialPerturbationIDs = [];
end
if nargin < 8
    event_func_handle = [];
end

% Initial STM is identity matrix (6x6)
STM0 = eye(6);
STM0_vec = reshape(STM0, 36, 1);

% Augmented initial state: [r; v; STM_vector]
x0_augmented = [rvec; vvec; STM0_vec];

% Augmented equations of motion
F = @(t, x) nBodyEOMWithSTM(t, x, muCB, mu_TBs, ephem_func, specialPerturbationIDs);

% ODE options
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Modified event function to work with augmented state
if ~isempty(event_func_handle)
    augmented_event = @(t, x) event_func_handle(t, x(1:6));
    options = odeset(options, 'Events', augmented_event);
    [tt, xx_augmented, te, ye_aug, ie] = ode45(F, timevector, x0_augmented, options);
    ye = ye_aug(:, 1:6); % Extract only position/velocity for events
else
    [tt, xx_augmented] = ode45(F, timevector, x0_augmented, options);
    te = []; ye = []; ie = [];
end

% Extract state and STM histories
yy = xx_augmented(:, 1:6);
STM_history = xx_augmented(:, 7:42);

end

function dxdt = nBodyEOMWithSTM(t, x, muCB, mu_TBs, ephem_func, specialPerturbationIDs)
    % Extract state and STM
    r_sc = x(1:3);
    v_sc = x(4:6);
    STM_vec = x(7:42);
    STM = reshape(STM_vec, 6, 6);
    
    % --- Calculate accelerations ---
    a_CB = -muCB * r_sc / norm(r_sc)^3;
    
    % Zonal harmonics
    a_zonal_pert = [0; 0; 0];
    if ~isempty(specialPerturbationIDs)
        j_term_types = abs(specialPerturbationIDs);
        all_zonal_accels = calculate_zonal_perturbations(r_sc, j_term_types);
        a_zonal_pert = sum(all_zonal_accels, 2);
    end
    
    % N-body perturbations with singularity avoidance
    a_nbody_pert = [0; 0; 0];
    min_distance_threshold1 = 10; % if set to 252, it does not find lower deltaV, no idea why
    
    if ~isempty(mu_TBs)
        r_TBs = ephem_func(t);
        for i = 1:length(mu_TBs)
            mu_tb = mu_TBs(i);
            r_tb = r_TBs(:, i);
            rtb2sc = r_sc - r_tb;
            dist_to_body = norm(rtb2sc);
            
            % Only include perturbation if we're not too close to the body
            % This prevents singularities when passing very close to Enceladus
            if dist_to_body > min_distance_threshold1
                pert_term_i = -mu_tb * (rtb2sc / dist_to_body^3) + mu_tb * (r_tb / norm(r_tb)^3);
                a_nbody_pert = a_nbody_pert + pert_term_i;
            end
            % If too close, the perturbation is effectively zero (we're in its SOI)
        end
    end
    
    a_total = a_CB + a_zonal_pert + a_nbody_pert;
    
    % --- Calculate A matrix (Jacobian of dynamics) ---
    A = computeAMatrix(r_sc, muCB, mu_TBs, ephem_func, t, specialPerturbationIDs);
    
    % --- STM differential equation: dΦ/dt = A * Φ ---
    STM_dot = A * STM;
    STM_dot_vec = reshape(STM_dot, 36, 1);
    
    % Assemble augmented state derivative
    dxdt = [v_sc; a_total; STM_dot_vec];
end

function A = computeAMatrix(r, muCB, mu_TBs, ephem_func, t, specialPerturbationIDs)
    % Compute the A matrix (Jacobian) for STM propagation
    % A = [∂ṙ/∂r  ∂ṙ/∂v]
    %     [∂v̇/∂r  ∂v̇/∂v]
    
    % Initialize A matrix
    A = zeros(6, 6);
    
    % Upper right block: ∂ṙ/∂v = I
    A(1:3, 4:6) = eye(3);
    
    % Lower right block: ∂v̇/∂v = 0 (for conservative forces)
    A(4:6, 4:6) = zeros(3, 3);
    
    % Lower left block: ∂v̇/∂r (gravity gradient matrix)
    
    % Central body contribution
    r_norm = norm(r);
    r_norm3 = r_norm^3;
    r_norm5 = r_norm^5;
    
    % Gravity gradient for point mass
    I3 = eye(3);
    rr_dyad = r * r';
    G_CB = -muCB/r_norm3 * I3 + 3*muCB/r_norm5 * rr_dyad;
    
    % Add J2 contribution if needed
    G_J2 = zeros(3, 3);
    if ~isempty(specialPerturbationIDs) && any(abs(specialPerturbationIDs) == 2)
        % G_J2 = computeJ2GravityGradient(r, muCB);
        J2_Saturn = 16290.573e-6;
        R_ref_Saturn = 60330;      % Equatorial radius for Saturn
        % G_J2 = computeJ2GravityGradient_generated(r, muCB, J2_Saturn, R_ref_Saturn);
        G_J2 = computeJ2GravityGradient_mathematical(r, muCB, J2_Saturn, R_ref_Saturn);
    end
    
    % N-body contributions with singularity avoidance
    G_nbody = zeros(3, 3);
    min_distance_threshold2 = 252; % do not crash into enceladus
    
    if ~isempty(mu_TBs)
        r_TBs = ephem_func(t);
        for i = 1:length(mu_TBs)
            mu_tb = mu_TBs(i);
            r_tb = r_TBs(:, i);
            r_rel = r - r_tb;
            r_rel_norm = norm(r_rel);
            
            % Only include gradient if not too close to the body
            if r_rel_norm > min_distance_threshold2
                r_rel_norm3 = r_rel_norm^3;
                r_rel_norm5 = r_rel_norm^5;
                
                G_tb = -mu_tb/r_rel_norm3 * I3 + 3*mu_tb/r_rel_norm5 * (r_rel * r_rel');
                G_nbody = G_nbody + G_tb;
            end
        end
    end
    
    A(4:6, 1:3) = G_CB + G_J2 + G_nbody;
end
