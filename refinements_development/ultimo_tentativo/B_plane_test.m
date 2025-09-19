clc;
clear all;
close all;

% this scripts checks that if there is no perturbers the nbody  propagator 
% propgates to the center of enceladus and its bplane paramters are all
% zero

%% ========================================================================
%  1. DEFINE INPUTS & INITIALIZE
%  ========================================================================

soi_multiplier = 64; 

pars.t_prop_hours  = 240; % hours
pars.GroundTr.npoints = 30e3;

pars.INPUTS.perturbingBodyNaifIDs = [10,602]; % no pertubers 
% pars.INPUTS.perturbingBodyNaifIDs = [-2, 602, 10]; % J2, Enceladus, Sun

% Perturbing bodies setup from the 'pars' struct
perturbingBodyNaifIDs = pars.INPUTS.perturbingBodyNaifIDs;
specialPerturbationIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs < 0);

pars.INPUTS.idCentral = 6;
pars.INPUTS.idMoon    = 1;
pars.INPUTS.Flyby.min_h = 10;

kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels);

% --- Load Gravitational Parameters ---
[pars.Planet.mu, ~, ~, ~] = planetConstants(pars.INPUTS.idCentral);
[pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = satMoonsConstants(pars.INPUTS.idMoon);
pars.Moon.HillSph = pars.Moon.OrbRad * (pars.Moon.mu / (3 * pars.Planet.mu))^(1/3);

% Dynamically build mu_TBs list
actualBodyNaifIDs = pars.INPUTS.perturbingBodyNaifIDs(pars.INPUTS.perturbingBodyNaifIDs >= 0);
mu_TBs = zeros(1, length(actualBodyNaifIDs));
for i = 1:length(actualBodyNaifIDs)
    id = actualBodyNaifIDs(i);
    switch id
        case 10,  mu_TBs(i) = getAstroConstants('Sun', 'Mu');
        case 5, mu_TBs(i) = getAstroConstants('Jupiter', 'Mu');
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

pars.INPUTS.NAIFCentral = 699;
pars.INPUTS.NAIFMoon = 602;

% propagate to 64 SOI
m_enceladus = 1.08e20;
m_saturn = getAstroConstants('Saturn', 'Mass');
a_enceladus = 238020;
r_soi_enceladus = a_enceladus * (m_enceladus / m_saturn)^(2/5);
pars.INPUTS.maxPropagationDistance = r_soi_enceladus * soi_multiplier;

spiceParam.frame    = 'J2000';
spiceParam.abcorr   = 'NONE';
spiceParam.observer = '699'; % Saturn Center NAIF ID


%%% flyby paramters
pars.INPUTS.epoch0 = date2mjd2000([2040 1 1 12 0 0]);
pars.INPUTS.V_inf = 4;
nodein = [pars.INPUTS.V_inf,0,0];
nodeout = [pars.INPUTS.V_inf,deg2rad(0.1), deg2rad(0.11)];



%% ========================================================================
%  2. CALCULATE LINKED-CONIC PERICENTER AND IN/OUT STATES
%  ========================================================================

% --- Get Enceladus' state at the flyby epoch ---
[r_enceladus_at_flyby, v_enceladus_at_flyby] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, pars.INPUTS.epoch0, true, spiceParam);

% --- Calculate incoming and outgoing V-infinity vectors ---
[vvinfin, r0_sc_in, v0_sc_in, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodein(1), nodein(2), nodein(3), pars.INPUTS.epoch0, pars);
[vvinfout, r0_sc_out, v0_sc_out, vvga] = vinfAlphaCrank_to_VinfCARTClaudio(nodeout(1), nodeout(2), nodeout(3), pars.INPUTS.epoch0, pars);

% --- Calculate pericenter state based on linked-conic assumptions ---
e_fly     = 1 + (( (pars.INPUTS.Flyby.min_h + pars.Moon.EquRad) * pars.INPUTS.V_inf^2) / pars.Moon.mu);
delta_max = 2 * asin(1 / e_fly);

Rm = buildRm(r_enceladus_at_flyby, v_enceladus_at_flyby);

[vvinfouBM, delta] = infVelBeforeDefect(vvinfin, vvinfout, delta_max, vvga);
vvinfin_bf = (Rm * vvinfin')';
vvinfouBM_bf = (Rm * vvinfouBM')';

Energy = 0.5 * norm(vvinfin_bf)^2;
sma = -pars.Moon.mu / (2 * Energy);
ecc = 1 / (sin(delta / 2));
rp = sma * (1 - ecc);
hhat = cross(vvinfin_bf, vvinfouBM_bf) / norm(cross(vvinfin_bf, vvinfouBM_bf));
vp = sqrt(norm(vvinfin_bf)^2 + 2 * pars.Moon.mu / rp);

rrp_bf = rp .* (vvinfin_bf - vvinfouBM_bf) / norm(vvinfin_bf - vvinfouBM_bf);
vvp_bf = vp .* cross(hhat, rrp_bf ./ rp);

Rm_inv = Rm';
rrp_saturn_centric = (Rm_inv * rrp_bf')' + r_enceladus_at_flyby;
vvp_saturn_centric = (Rm_inv * vvp_bf')' + v_enceladus_at_flyby;

pericenter_LC = [rrp_saturn_centric(:); vvp_saturn_centric(:)];

%% ========================================================================
%  3. PROPAGATE LINKED CONICS TO - 64 ENCELADUS SOI
%  ========================================================================

% --- Define time vectors for propagation ---
duration_sec = pars.t_prop_hours * 3600;
time_steps = pars.GroundTr.npoints;
time_vector_fwd = linspace(eps, duration_sec, time_steps)';
time_vector_bwd = linspace(0, -duration_sec, time_steps)';

% --- Setup event function to stop at 64 SOI distance ---
soi_64_dist = pars.INPUTS.maxPropagationDistance;
ephem_enceladus = @(t) get_body_positions_wrapper(t, pars.INPUTS.epoch0, pars.INPUTS.NAIFMoon, spiceParam);
soi_event_func = @(t, x) soiCrossingEvent(t, x, ephem_enceladus, soi_64_dist);

% --- Propagate LC backward the -64 SOI boundary ---
[timeLC_bwd, stateLC_bwd, te_bwd, ye_bwd, ~] = propagateKeplerODE2(r0_sc_in', v0_sc_in', time_vector_bwd, pars.Planet.mu, soi_event_func);

% --- Propagate LC forward the +64 SOI boundary --
[timeLC_fwd, stateLC_fwd, te_fwd, ye_fwd, ~] = propagateKeplerODE2(r0_sc_out', v0_sc_out', time_vector_fwd,  pars.Planet.mu, soi_event_func);


LC_in = ye_bwd(1, :)';

%% ========================================================================
%  4. PROPAGATE N-BODY FROM LC_in
%  ========================================================================

% --- Set initial state for N-body propagation ---
NB_in = LC_in;
initial_epoch_offset_days = pars.INPUTS.epoch0 + te_bwd/86400; % Time offset from pericenter epoch

% --- Propagate forward from NB_in ---
flyby_tot_real_duration_seconds =  (-te_bwd+te_fwd); 
time_vector_nbody = linspace(eps, flyby_tot_real_duration_seconds*1.5, time_steps)'; % Adjusted time vector

% Define ephemeris handle for the integrator
ephem_handle = @(t_sec) get_body_positions_wrapper(t_sec, initial_epoch_offset_days, actualBodyNaifIDs, spiceParam);

% Redefine event function handle for the N-body start time
ephem_enceladus_nbody = @(t) get_body_positions_wrapper(t, initial_epoch_offset_days, pars.INPUTS.NAIFMoon, spiceParam);

% soi event stops when arriving at +64 SOI(fgoing away from Enceladu in the
% other direction)
soi_event_func_nbody = @(t, x) soiExitEvent(t, x, ephem_enceladus_nbody, soi_64_dist); % <-- CHANGE HERE

[time_out_nb, state_out_nb, te_fwd_nb, ye_fwd_nb, ~] = propagateNBodyODE3(NB_in(1:3), NB_in(4:6), ...
    time_vector_nbody, pars.Planet.mu, mu_TBs, ephem_handle, specialPerturbationIDs, soi_event_func_nbody);

%% ========================================================================
%  5. FIND N-BODY PERICENTER
%  ========================================================================

% --- Calculate distance to Enceladus for the entire N-body trajectory ---
num_steps_nb = length(time_out_nb);
dist_to_enceladus = zeros(num_steps_nb, 1);
for i = 1:num_steps_nb
    % Get Enceladus position at the corresponding time
    t_sec_from_start_of_nbody = time_out_nb(i);
    r_enc_current = get_body_positions_wrapper(t_sec_from_start_of_nbody, initial_epoch_offset_days, pars.INPUTS.NAIFMoon, spiceParam);
    
    % Calculate distance
    dist_to_enceladus(i) = norm(state_out_nb(i, 1:3) - r_enc_current');
end

% --- Find the minimum distance to locate the pericenter ---
[~, idx_pericenter] = min(dist_to_enceladus);
min_dist_nb = min(dist_to_enceladus);

pericenter_NB = state_out_nb(idx_pericenter, :)';
time_at_pericenter_sec_from_NB_start = time_out_nb(idx_pericenter);
epoch_pericenter_NB = initial_epoch_offset_days + time_at_pericenter_sec_from_NB_start / 86400;



%% ========================================================================
%  6. PRINT THE RESULTS
%  ========================================================================

% --- Confront N-Body and Linked-Conic Pericenters ---
fprintf('Pericenter State Comparison (Saturn-Centric J2000):\n');
fprintf('------------------------------------------------------\n');
fprintf('Linked Conic Pericenter State:\n');
fprintf('r = [%16.6f, %16.6f, %16.6f] km\n', pericenter_LC(1:3));
fprintf('v = [%16.6f, %16.6f, %16.6f] km/s\n\n', pericenter_LC(4:6));

fprintf('N-Body Pericenter State:\n');
fprintf('r = [%16.6f, %16.6f, %16.6f] km\n', pericenter_NB(1:3));
fprintf('v = [%16.6f, %16.6f, %16.6f] km/s\n\n', pericenter_NB(4:6));

% Calculate differences
pos_diff = pericenter_NB(1:3) - pericenter_LC(1:3);
vel_diff = pericenter_NB(4:6) - pericenter_LC(4:6);
altitude_LC = rp - pars.Moon.EquRad;
altitude_NB = min_dist_nb - pars.Moon.EquRad;

fprintf('Difference Analysis:\n');
fprintf('--------------------\n');
fprintf('Position difference norm: %.6f km\n', norm(pos_diff));
fprintf('Velocity difference norm: %.6f km/s\n\n', norm(vel_diff));

fprintf('Pericenter Altitude Comparison (relative to Enceladus):\n');
fprintf('-----------------------------------------------------\n');
fprintf('Linked Conic Altitude: %.6f km\n', altitude_LC);
fprintf('N-Body Altitude:       %.6f km\n', altitude_NB);
fprintf('Altitude Difference:   %.6f km\n', altitude_NB - altitude_LC);
fprintf('Enceladus Radius:   %.6f km\n', pars.Moon.EquRad);
fprintf('========================================================================\n\n');

%% ========================================================================
%  7. B-PLANE ANALYSIS
%  ========================================================================

fprintf('B-Plane Targeting Analysis:\n');
fprintf('---------------------------\n');

% --- Step A: Calculate the TARGET B-Plane from Linked-Conics Asymptotes ---

% Use the original V-infinity vectors that define the desired 3D geometry
[B_vec_target, rp_target_from_asymptotes, B_R_target, B_T_target] = ...
    b_plane_from_asymptotes(vvinfin', vvinfouBM', pars.Moon.mu);

fprintf('Target B-Plane Parameters (from Linked Conics):\n');
fprintf('  - B_R Target: %.6f km\n', B_R_target);
fprintf('  - B_T Target: %.6f km\n', B_T_target);
fprintf('  - |B| Target: %.6f km\n\n', norm(B_vec_target));


% --- Step B: Calculate the ACHIEVED B-Plane from N-Body Pericenter ---

% First, get Enceladus's state at the moment of N-body pericenter
[r_enc_at_nb_peri, v_enc_at_nb_peri] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, epoch_pericenter_NB, true, spiceParam);

% Calculate the spacecraft state relative to Enceladus
r_sc_rel_enc = pericenter_NB(1:3) - r_enc_at_nb_peri';
v_sc_rel_enc = pericenter_NB(4:6) - v_enc_at_nb_peri';

% Calculate the square of the hyperbolic excess velocity (v_inf^2)
% from the energy of the relative orbit: v_inf^2 = v^2 - 2*mu/r
v_inf_sq_achieved = dot(v_sc_rel_enc, v_sc_rel_enc) - 2 * pars.Moon.mu / norm(r_sc_rel_enc);

% Now, calculate the achieved B-plane parameters
[B_R_achieved, B_T_achieved, B_vec_achieved] = ...
    b_plane_targeting(r_sc_rel_enc, v_sc_rel_enc, v_inf_sq_achieved, pars.Moon.mu);

fprintf('Achieved B-Plane Parameters (from N-Body propagation):\n');
fprintf('  - B_R Achieved: %.6f km\n', B_R_achieved);
fprintf('  - B_T Achieved: %.6f km\n', B_T_achieved);
fprintf('  - |B| Achieved: %.6f km\n\n', norm(B_vec_achieved));

% --- Step C: Calculate the B-Plane Difference ("Miss Vector") ---

% The miss vector is the difference between the achieved and target B-vectors
B_miss_vec = B_vec_achieved - B_vec_target;
miss_distance = norm(B_miss_vec);

% Differences in the components
diff_B_R = B_R_achieved - B_R_target;
diff_B_T = B_T_achieved - B_T_target;

fprintf('B-Plane Miss Analysis:\n');
fprintf('  - Miss Distance (|B_achieved - B_target|): %.6f km\n', miss_distance);
fprintf('  - Delta B_R (Error in R-component):      %.6f km\n', diff_B_R);
fprintf('  - Delta B_T (Error in T-component):      %.6f km\n', diff_B_T);
fprintf('========================================================================\n\n');


% --- Verification Step: Confirming Geometric Consistency ---
fprintf('Verification of Geometric Consistency:\n');
fprintf('  - Linked-Conics Model rp: %.6f km\n', rp);
fprintf('  - B-Plane Target rp (from asymptotes): %.6f km\n\n', rp_target_from_asymptotes);


%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function positions = get_body_positions_wrapper(t_sec, t0_mjd_start, naif_ids, sp_param)
    current_mjd = t0_mjd_start + t_sec / 86400;
    num_bodies = length(naif_ids);
    positions = zeros(3, num_bodies);

    for k = 1:num_bodies
        [r, ~] = EphSS_car_spice2(naif_ids(k), current_mjd, true, sp_param);
        positions(:, k) = r(:); % Ensure column vector
    end
end


function [value, isterminal, direction] = soiCrossingEvent(t, x, ephem_enceladus_handle, target_distance)
    r_sc = x(1:3);
    r_enceladus = ephem_enceladus_handle(t);
    current_distance = norm(r_sc - r_enceladus);
    value = current_distance - target_distance;
    isterminal = 1;
    direction = 0;
end


function [value, isterminal, direction] = soiExitEvent(t, x, ephem_enceladus_handle, target_distance)
    % Helper function to detect when the spacecraft exits a target distance boundary.
    
    r_sc = x(1:3); % Spacecraft position vector
    r_enceladus = ephem_enceladus_handle(t); % Enceladus position vector
    
    current_distance = norm(r_sc - r_enceladus);
    
    % The event occurs when the current distance equals the target distance.
    value = current_distance - target_distance;
    
    % Terminate the integration when the event occurs.
    isterminal = 1;
    
    % Trigger only when the event function is increasing (i.e., distance is increasing).
    direction = 1;
end