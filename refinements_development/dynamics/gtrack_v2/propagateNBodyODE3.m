function [tt, yy, te, ye, ie] = propagateNBodyODE3(rvec, vvec, timevector, muCB, mu_TBs, ephem_func, specialPerturbationIDs, event_func_handle)
%
%   Inputs
%   ======
%   rvec        3x1 initial position of spacecraft w.r.t. central body (km)
%   vvec        3x1 initial velocity of spacecraft w.r.t. central body (km/s)
%   timevector  Nx1 vector of time points (s)
%   muCB        Scalar gravitational parameter of central body (km^3/s^2)
%   mu_TBs      Mx1 vector of gravitational parameters for M target bodies
%   ephem_func  Function handle that takes time 't' and returns a 3xM
%               matrix of target body positions w.r.t. the central body.
%               Example: r_moons = ephem_func(t);
%   event_func_handle   (handle) Optional function handle for event detection.
%
%   Outputs
%   =======
%    tt          Nx1 output times (s)
%   yy           Nx6 matrix of state history [x, y, z, vx, vy, vz]
%   te, ye, ie   (vector) Event time, state, and index.

if nargin < 7
    specialPerturbationIDs = []; % Default to no special perturbations
end
if nargin < 8
    event_func_handle = []; % Default to no events
end

% --- Pass all parameters to the Equations of Motion (EOM) ---
F = @(t, x) nBodyEOM(t, x, muCB, mu_TBs, ephem_func, specialPerturbationIDs);

% --- event function ---
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
if ~isempty(event_func_handle)
    options = odeset(options, 'Events', event_func_handle);
end

   %  check if an event function was provided.
    % If it was, we add it to the options and call ode45 requesting all 5 outputs.
    % If not, we call ode45 requesting only 2 outputs and set the event
    % outputs to empty arrays to maintain a consistent function signature.
if ~isempty(event_func_handle)
        options = odeset(options, 'Events', event_func_handle);
        [tt, yy, te, ye, ie] = ode45(F, timevector, [rvec; vvec], options);
    else
        [tt, yy] = ode45(F, timevector, [rvec; vvec], options);
        % Define empty event outputs when no event function is used
        te = [];
        ye = [];
        ie = [];
end

% --- Nested function for the Equations of Motion (unchanged) ---
function dxdt = nBodyEOM(t, x, muCB, mu_TBs, ephem_func, specialPerturbationIDs)
    r_sc = x(1:3); % Spacecraft position vector

    % 1. Acceleration from the Central Body (point mass)
    a_CB = -muCB * r_sc / (norm(r_sc)^3);

    % 2. Acceleration from Zonal Harmonic Perturbations (J2, J3, etc.)
    a_zonal_pert = [0; 0; 0];
    if ~isempty(specialPerturbationIDs)
        j_term_types = abs(specialPerturbationIDs);
        all_zonal_accels = calculate_zonal_perturbations(r_sc, j_term_types);
        a_zonal_pert = sum(all_zonal_accels, 2);
    end

    % 3. Acceleration from N-body perturbers
    a_nbody_pert = [0; 0; 0];
    if ~isempty(mu_TBs)
        r_TBs = ephem_func(t); % Get positions of all perturbing bodies
        for i = 1:length(mu_TBs)
            mu_tb = mu_TBs(i);
            r_tb  = r_TBs(:, i);
            rtb2sc = r_sc - r_tb;
            dist_to_body = norm(rtb2sc);
            
            if dist_to_body > 10
                pert_term_i = -mu_tb * (rtb2sc / dist_to_body^3) + mu_tb * (r_tb / norm(r_tb)^3);
                a_nbody_pert = a_nbody_pert + pert_term_i;
            end
        end
    end
    
    % 4. Total Acceleration is the sum of all forces
    a_total = a_CB + a_zonal_pert + a_nbody_pert;
    
    % 5. Assemble the state derivative vector
    dxdt = [x(4:6); a_total];
end

end