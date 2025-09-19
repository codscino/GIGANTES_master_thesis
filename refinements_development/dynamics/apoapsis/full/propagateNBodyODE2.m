function [tt, yy, te, ye, ie] = propagateNBodyODE2(rvec, vvec, timevector, muCB, mu_TBs, ephem_func, specialPerturbationIDs, event_func)
%
% Propagate a point-mass in an n-body system with optional event detection.
%
% Inputs
% ======
% rvec 3x1 initial position of spacecraft w.r.t. central body (km)
% vvec 3x1 initial velocity of spacecraft w.r.t. central body (km/s)
% timevector Nx1 vector of time points (s)
% muCB Scalar gravitational parameter of central body (km^3/s^2)
% mu_TBs Mx1 vector of gravitational parameters for M target bodies
% ephem_func Function handle that takes time 't' and returns a 3xM
%            matrix of target body positions w.r.t. the central body.
%            Example: r_moons = ephem_func(t);
% specialPerturbationIDs (vector) List of special perturbation IDs (e.g., -2 for J2).
% event_func (optional) Event function handle for event detection
%            Should return [value, isterminal, direction]
%
% Outputs
% =======
% tt Time vector
% yy State vector history
% te Event times (empty if no event function provided)
% ye Event states (empty if no event function provided)
% ie Event indices (empty if no event function provided)

% Handle optional parameters
if nargin < 7
    specialPerturbationIDs = []; % Default to no special perturbations
end

if nargin < 8
    event_func = []; % Default to no event function
end

% --- Pass all parameters to the Equations of Motion (EOM) ---
F = @(t, x) nBodyEOM(t, x, muCB, mu_TBs, ephem_func, specialPerturbationIDs);

% Set up options
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Add event function if provided
if ~isempty(event_func)
    options = odeset(options, 'Events', event_func);
end

% Run integration
if ~isempty(event_func)
    [tt, yy, te, ye, ie] = ode45(F, timevector, [rvec; vvec], options);
else
    [tt, yy] = ode45(F, timevector, [rvec; vvec], options);
    te = [];
    ye = [];
    ie = [];
end

% --- Nested function for the Equations of Motion ---
function dxdt = nBodyEOM(t, x, muCB, mu_TBs, ephem_func, specialPerturbationIDs)
    r_sc = x(1:3); % Spacecraft position vector

    % 1. Acceleration from the Central Body (point mass)
    a_CB = -muCB * r_sc / (norm(r_sc)^3);

    % 2. Acceleration from Zonal Harmonic Perturbations (J2, J3, etc.)
    a_zonal_pert = [0; 0; 0];
    if ~isempty(specialPerturbationIDs)
        % Convert negative special IDs to positive J-term numbers (e.g., [-2, -4] -> [2, 4])
        j_term_types = abs(specialPerturbationIDs);
        % Call the zonal calculation function ONCE to get all requested accelerations
        all_zonal_accels = calculate_zonal_perturbations(r_sc, j_term_types);
        % Sum the columns of the output matrix to get the total zonal perturbation
        a_zonal_pert = sum(all_zonal_accels, 2);
    end

    % 3. Acceleration from N-body perturbers
    a_nbody_pert = [0; 0; 0];
    if ~isempty(mu_TBs)
        r_TBs = ephem_func(t); % Get positions of all perturbing bodies
        for i = 1:length(mu_TBs)
            mu_tb = mu_TBs(i);
            r_tb = r_TBs(:, i);
            rtb2sc = r_sc - r_tb;
            pert_term_i = -mu_tb * (rtb2sc / norm(rtb2sc)^3) + mu_tb * (r_tb / norm(r_tb)^3);
            a_nbody_pert = a_nbody_pert + pert_term_i;
        end
    end

    % 4. Total Acceleration is the sum of all forces
    a_total = a_CB + a_zonal_pert + a_nbody_pert;

    % 5. Assemble the state derivative vector
    dxdt = [x(4:6); a_total];
end

end