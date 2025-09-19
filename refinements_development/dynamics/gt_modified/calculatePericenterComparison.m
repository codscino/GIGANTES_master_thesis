function [rp_target, rp_nbody] = calculatePericenterComparison(pars, spiceParam, t0_mjd, perturbingBodyNaifIDs, useJ2)
% Compares a targeted flyby pericenter with
% the result from a high-fidelity N-body propagation.

if nargin < 5
    useJ2 = false;
end

%% ========================================================================
%  1. DEFINE TARGET PERICENTER USING EXISTING LINKED CONICS FUNCTIONS
%  ========================================================================

% --- Load necessary constants ---
[pars.Planet.mu, ~, ~, ~] = planetConstants(pars.INPUTS.idCentral);
[pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = satMoonsConstants(pars.INPUTS.idMoon);

% --- Define flyby nodes from input parameters ---
% The design of the pericenter depends on the incoming and outgoing asymptotes.
nodein  = [pars.INPUTS.V_inf, pars.INPUTS.alpha, pars.INPUTS.k];

% dummy outgoing crank angle -> is this ok ???
kou = pars.INPUTS.k + deg2rad(1); 
nodeout = [pars.INPUTS.V_inf, pars.INPUTS.alpha, kou];

% --- Get inertial V-infinity vectors and moon state at the target time ---
[vvinfin, ~, ~, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodein(1), nodein(2), nodein(3), t0_mjd);
[vvinfou, rr_moon_t0, vv_moon_t0, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodeout(1), nodeout(2), nodeout(3), t0_mjd);

% --- Calculate delta_max ---
rp_target_mag = pars.Moon.EquRad + pars.Moon.hmin;
e_fly     = 1 + ((rp_target_mag * pars.INPUTS.V_inf^2)/pars.Moon.mu);
pars.delta_max = 2 * asin(1 / e_fly);

% --- get the ZSOI/Target pericenter info ---
[~, ~, ~, ~, ~, ~, ~, ~, lat_deg_target, lon_deg_target, rrp_moon_fixed, vvp_moon_fixed] = ...
    groundTrackFlyByMoonClaudio(vvinfin, vvinfou, pars.delta_max, rr_moon_t0, vv_moon_t0, pars.Moon.mu, pars);

% --- Store Target Pericenter Info ---
rp_target.lat_deg = lat_deg_target;
rp_target.lon_deg = lon_deg_target;
rp_target.alt_km  = norm(rrp_moon_fixed) - pars.Moon.EquRad;

% --- Convert pericenter state from body fixed to ICRF 
R_enc = 252.3; % average from spice
r_enc = EphSS_car_spice2(602, t0_mjd, true);
rrp_moon_icrf = enceladus2icrf(rrp_moon_fixed, r_enc, R_enc, t0_mjd);
vvp_moon_icrf = enceladus2icrf(vvp_moon_fixed, r_enc, R_enc, t0_mjd);

[r_sat, v_sat, ~] = EphSS_car_spice2(699, t0_mjd, true);
% convert to column vectors
rrp_moon_icrf = (rrp_moon_icrf-r_sat)';
vvp_moon_icrf = (vvp_moon_icrf- v_sat)';

% b1 = -rr_moon_t0./norm(rr_moon_t0);
% b3 = cross(rr_moon_t0, vv_moon_t0)./norm(cross(rr_moon_t0, vv_moon_t0));
% b2 = cross(b3,b1);
% % This is the matrix that rotates from ICRF to the moon-fixed frame.
% Rm_I2M = [b1; b2; b3];
% % The transpose rotates from moon-fixed back to ICRF.
% rrp_moon_icrf = Rm_I2M' * rrp_moon_fixed';
% vvp_moon_icrf = Rm_I2M' * vvp_moon_fixed';


%% ========================================================================
%  2. BACK-PROPAGATE TO FIND THE COMMON STARTING POINT
%  ========================================================================

% --- Calculate SOI and time to propagate ---
soi_radius = pars.Moon.OrbRad * (pars.Moon.mu / (3 * pars.Planet.mu))^(1/3);
a_hyperbola = -pars.Moon.mu / (pars.INPUTS.V_inf^2);
ecc_hyperbola = 1 + (norm(rrp_moon_fixed) * pars.INPUTS.V_inf^2) / pars.Moon.mu;

% Eccentric Anomaly (H) at SOI boundary
cosh_H_soi = (soi_radius - a_hyperbola) / (a_hyperbola * ecc_hyperbola);
H_soi = acosh(max(1, cosh_H_soi));

% Time from pericenter to SOI using Kepler's Equation for hyperbolas
n_hyperbola = sqrt(-pars.Moon.mu / a_hyperbola^3);
dt_soi = (1/n_hyperbola) * (ecc_hyperbola * sinh(H_soi) - H_soi);

% --- Propagate backward from target pericenter to SOI edge ---
[r_soi_moon, v_soi_moon] = propagateKeplerUV(rrp_moon_icrf, vvp_moon_icrf, -dt_soi, pars.Moon.mu);

% --- Define the common starting point in the planet-centric inertial frame ---
entry_time_mjd = t0_mjd - dt_soi / 86400;
[rr_moon_entry, vv_moon_entry, ~] = EphSS_car_spice2(602, entry_time_mjd, true, spiceParam);
%[rr_moon_entry, vv_moon_entry, ~] = approxEphem_CC(pars.INPUTS.idMoon, entry_time_mjd, pars.INPUTS.idCentral);

r0_common = rr_moon_entry + r_soi_moon;
v0_common = vv_moon_entry + v_soi_moon;


%% ========================================================================
%  3. FORWARD-PROPAGATE N-BODY AND FIND ACTUAL PERICENTER
%  ========================================================================

% --- Set up N-body propagation ---
simulation_duration_sec = 2.5 * dt_soi; % Propagate well past the exit
time_vector = linspace(0, simulation_duration_sec, 2000)'; % High resolution for accuracy

% --- Dynamically build the list of gravitational parameters (mu) ---
mu_TBs = zeros(1, length(perturbingBodyNaifIDs));
for i = 1:length(perturbingBodyNaifIDs)
    id = perturbingBodyNaifIDs(i);
    switch id
        case 10 % Sun
            mu_TBs(i) = getAstroConstants('Sun', 'Mu');
        case 599 % Jupiter
            mu_TBs(i) = getAstroConstants('Jupiter', 'Mu');
        case 601 % Mimas
            [~, mu, ~, ~] = satMoonsConstants(0); mu_TBs(i) = mu;
        case 602 % Enceladus
            [~, mu, ~, ~] = satMoonsConstants(1); mu_TBs(i) = mu;
        case 603 % Tethys
            [~, mu, ~, ~] = satMoonsConstants(2); mu_TBs(i) = mu;
        case 604 % Dione
            [~, mu, ~, ~] = satMoonsConstants(3); mu_TBs(i) = mu;
        case 605 % Rhea
            [~, mu, ~, ~] = satMoonsConstants(4); mu_TBs(i) = mu;
        case 606 % Titan
            [~, mu, ~, ~] = satMoonsConstants(5); mu_TBs(i) = mu;
        case 607 % Hyperion
            mu_TBs(i) = 0.374; % From script
        case 608 % Iapetus
            mu_TBs(i) = 120.4; % From script
        otherwise
            error('calculateApoapsisComparison:UnknownNaifID', ...
                  'The NAIF ID %d is not recognized for mu calculation.', id);
    end
end

ephem_handle = @(t_sec) get_all_moon_positions(t_sec, entry_time_mjd, perturbingBodyNaifIDs, spiceParam);

% --- Run the N-body simulation ---
[~, nbody_state] = propagateNBodyODE(r0_common', v0_common', time_vector, ...
                                   pars.Planet.mu, mu_TBs, ephem_handle, useJ2);

% --- Find the point of closest approach in the N-body simulation ---
sim_times_mjd = entry_time_mjd + time_vector/86400;
moon_pos_history = zeros(length(time_vector), 3);
for i = 1:length(time_vector)
    [moon_pos_history(i,:), ~, ~] = EphSS_car_spice2(602, sim_times_mjd(i), true, spiceParam);
    %[moon_pos_history(i,:), ~, ~] = approxEphem_CC(pars.INPUTS.idMoon, sim_times_mjd(i), pars.INPUTS.idCentral);
end
relative_dist_nbody = vecnorm(nbody_state(:, 1:3) - moon_pos_history, 2, 2);
[min_dist_nbody, idx_min_nbody] = min(relative_dist_nbody);

% --- Store N-Body Pericenter Info ---
rp_nbody.alt_km = min_dist_nbody - pars.Moon.EquRad;
time_at_nbody_rp_mjd = sim_times_mjd(idx_min_nbody);

% Position of spacecraft and moon at the moment of N-body closest approach
r_sc_at_rp_icrf = nbody_state(idx_min_nbody, 1:3);
r_moon_at_rp_icrf = moon_pos_history(idx_min_nbody, :);

% Convert the N-body pericenter position to lat/lon
r_iau_nbody = icrf2enceladus(r_sc_at_rp_icrf, r_moon_at_rp_icrf, time_at_nbody_rp_mjd);
rp_nbody.lat_deg = asind(r_iau_nbody(3));
rp_nbody.lon_deg = atan2d(r_iau_nbody(2), r_iau_nbody(1));

end


% =========================================================================
%                       HELPER FUNCTIONS BELOW
% =========================================================================

function [rf, vf] = propagateKeplerUV(r0, v0, dt, mu)
% Propagates a state vector using a universal variable formulation,
% robust for elliptical, parabolic, and hyperbolic orbits.

    r0_norm = norm(r0);
    v0_norm = norm(v0);

    % Orbital parameters
    vr0 = dot(r0, v0) / r0_norm;
    alpha = 2/r0_norm - v0_norm^2/mu; % Reciprocal of semi-major axis (1/a)

    % Initial guess for universal variable X
    X = sqrt(mu) * abs(alpha) * dt;

    % Newton-Raphson iteration to solve Stumpff's universal Kepler equation
    n_iter = 0;
    max_iter = 10;
    tol = 1e-8;
    ratio = 1;

    while abs(ratio) > tol && n_iter < max_iter
        z = alpha * X^2;
        [C, S] = stumpff(z);
        
        f = (r0_norm * vr0 / sqrt(mu)) * X^2 * C + (1 - alpha * r0_norm) * X^3 * S + r0_norm * X - sqrt(mu) * dt;
        df_dX = (r0_norm * vr0 / sqrt(mu)) * X * (1 - alpha * X^2 * S) + (1 - alpha * r0_norm) * X^2 * C + r0_norm;
        
        ratio = f / df_dX;
        X = X - ratio;
        n_iter = n_iter + 1;
    end

    % Lagrange F and G coefficients
    z = alpha * X^2;
    [C, S] = stumpff(z);
    
    f_lagrange = 1 - (X^2 / r0_norm) * C;
    g_lagrange = dt - (1/sqrt(mu)) * X^3 * S;

    % Final state vector
    rf = (f_lagrange * r0 + g_lagrange * v0)';
    
    rf_norm = norm(rf);
    df_dt_lagrange = (sqrt(mu) / (rf_norm * r0_norm)) * (alpha * X^3 * S - X);
    g_dot_lagrange = 1 - (X^2 / rf_norm) * C;
    
    vf = (df_dt_lagrange * r0 + g_dot_lagrange * v0)';
end

function [C, S] = stumpff(z)
% Computes Stumpff functions C(z) and S(z).
    if z > 1e-6
        S = (sqrt(z) - sin(sqrt(z))) / (sqrt(z))^3;
        C = (1 - cos(sqrt(z))) / z;
    elseif z < -1e-6
        S = (sinh(sqrt(-z)) - sqrt(-z)) / (sqrt(-z))^3;
        C = (cosh(sqrt(-z)) - 1) / (-z);
    else
        S = 1/6 - z/120 + z^2/5040;
        C = 1/2 - z/24 + z^2/720;
    end
end
