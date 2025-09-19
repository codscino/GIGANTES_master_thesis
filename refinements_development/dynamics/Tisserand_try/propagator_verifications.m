% clc;
clear all;
close all;

%% ========================================================================
%  1. DEFINE INPUTS & INITIALIZE
%  ========================================================================

% soi_multiplier = 136.6; 
soi_multiplier = 64; 

pars.GroundTr.t_prop  = 50*2000;
pars.GroundTr.npoints = 30e3;

pars.INPUTS.perturbingBodyNaifIDs = [-2,10,602];
% pars.INPUTS.perturbingBodyNaifIDs = [];


pars.INPUTS.idCentral = 6;
pars.INPUTS.idMoon    = 1;
pars.INPUTS.Flyby.min_h = 25;

kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels);

% --- Load Gravitational Parameters ---
[pars.Planet.mu, ~, ~, ~] = planetConstants(pars.INPUTS.idCentral);
[pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = satMoonsConstants(pars.INPUTS.idMoon);
pars.Moon.HillSph = pars.Moon.OrbRad * (pars.Moon.mu / (3 * pars.Planet.mu))^(1/3);

% Dynamically build mu_TBs list
actualBodyNaifIDs_init = pars.INPUTS.perturbingBodyNaifIDs(pars.INPUTS.perturbingBodyNaifIDs >= 0);
mu_TBs = zeros(1, length(actualBodyNaifIDs_init));
for i = 1:length(actualBodyNaifIDs_init)
    id = actualBodyNaifIDs_init(i);
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
pars.INPUTS.mu_TBs = mu_TBs;

pars.INPUTS.NAIFCentral = 699;
pars.INPUTS.NAIFMoon = 602;

% propagate to 64 SOI
m_enceladus = 1.08e20;
m_saturn = getAstroConstants('Saturn', 'Mass');
a_enceladus = 238020;
r_soi_enceladus = a_enceladus * (m_enceladus / m_saturn)^(2/5);
pars.INPUTS.maxPropagationDistance = r_soi_enceladus * soi_multiplier;

% --- Set baseline constant parameters for loops ---
BASELINE_VINF = 4;

pars.INPUTS.V_inf = BASELINE_VINF;

%% ========================================================================
%  2A. EPOCH LOOP (MUST RUN FIRST TO FIND OPTIMAL EPOCH)
%  ========================================================================

% Use a non-zero deflection for a realistic scenario
nodein_base = [BASELINE_VINF, 0, 0];
% nodeout_base = [BASELINE_VINF, 0.15, deg2rad(1)]; 
nodeout_base = [BASELINE_VINF, 0.15, deg2rad(1)]; 

pars.INPUTS.epoch0 = date2mjd2000([2040, 1, 1, 0, 0, 0]);

%% ====================================================================
%  START: Dismantled propagateAndCompare (Sections 1-3)
%  ====================================================================

%% Section 1: INITIALIZE PARAMETERS & SPICE
%  --------------------------------------------------------------------
spiceParam.frame    = 'J2000';
spiceParam.abcorr   = 'NONE';
spiceParam.observer = '699'; % Saturn Center NAIF ID

mu_central_body = pars.Planet.mu;
mu_enceladus = pars.Moon.mu;

% Perturbing bodies setup from the 'pars' struct
perturbingBodyNaifIDs = pars.INPUTS.perturbingBodyNaifIDs;
actualBodyNaifIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs >= 0);
specialPerturbationIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs < 0);
mu_TBs_local = pars.INPUTS.mu_TBs;

%% Section 2: CALCULATE LINKED-CONIC PERICENTER AND IN/OUT STATES
%  --------------------------------------------------------------------

% --- Get Enceladus' state at the flyby epoch ---
[r_enceladus_at_flyby, v_enceladus_at_flyby] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, pars.INPUTS.epoch0, true, spiceParam);

% --- Calculate incoming and outgoing V-infinity vectors ---
[vvinfin, r0_sc_in, v0_sc_in, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodein_base(1), nodein_base(2), nodein_base(3), pars.INPUTS.epoch0, pars);
[vvinfout, r0_sc_out, v0_sc_out, vvga] = vinfAlphaCrank_to_VinfCARTClaudio(nodeout_base(1), nodeout_base(2), nodeout_base(3), pars.INPUTS.epoch0, pars);

% --- Calculate pericenter state based on linked-conic assumptions ---
e_fly     = 1 + (( (pars.INPUTS.Flyby.min_h + pars.Moon.EquRad) * pars.INPUTS.V_inf^2) / mu_enceladus);
delta_max = 2 * asin(1 / e_fly);

Rm = buildRm(r_enceladus_at_flyby, v_enceladus_at_flyby);

[vvinfouBM, delta] = infVelBeforeDefect(vvinfin, vvinfout, delta_max, vvga);
vvinfin_bf = (Rm * vvinfin')';
vvinfouBM_bf = (Rm * vvinfouBM')';

Energy = 0.5 * norm(vvinfin_bf)^2;
sma = -mu_enceladus / (2 * Energy);
ecc = 1 / (sin(delta / 2));
rp = sma * (1 - ecc);
hhat = cross(vvinfin_bf, vvinfouBM_bf) / norm(cross(vvinfin_bf, vvinfouBM_bf));
vp = sqrt(norm(vvinfin_bf)^2 + 2 * mu_enceladus / rp);

rrp_bf = rp .* (vvinfin_bf - vvinfouBM_bf) / norm(vvinfin_bf - vvinfouBM_bf);
vvp_bf = vp .* cross(hhat, rrp_bf ./ rp);

Rm_inv = Rm';
rrp_saturn_centric = (Rm_inv * rrp_bf')' + r_enceladus_at_flyby;
vvp_saturn_centric = (Rm_inv * vvp_bf')' + v_enceladus_at_flyby;

pericenter_LC = [rrp_saturn_centric(:); vvp_saturn_centric(:)];

%% Section 3: PROPAGATE LINKED CONICS TO +/- 64 ENCELADUS SOI
%  --------------------------------------------------------------------

% --- Define time vectors for propagation ---
duration_sec = pars.GroundTr.t_prop * 60;
time_steps = pars.GroundTr.npoints;
time_vector_fwd = linspace(eps, duration_sec, time_steps)';
time_vector_bwd = linspace(0, -duration_sec, time_steps)';

% --- Setup event function to stop at 64 SOI distance ---
soi_64_dist = pars.INPUTS.maxPropagationDistance;
ephem_enceladus = @(t) get_body_positions_wrapper(t, pars.INPUTS.epoch0, pars.INPUTS.NAIFMoon, spiceParam);
soi_event_func = @(t, x) soiCrossingEvent(t, x, ephem_enceladus, soi_64_dist);

% --- Propagate LC backward and forward to the 64 SOI boundary ---
[~, ~, te_bwd, ye_bwd, ~] = propagateKeplerODE2(r0_sc_in', v0_sc_in', time_vector_bwd, mu_central_body, soi_event_func);
[~, ~, te_fwd, ye_fwd, ~] = propagateKeplerODE2(r0_sc_out', v0_sc_out', time_vector_fwd, mu_central_body, soi_event_func);

if isempty(ye_fwd) || isempty(ye_bwd)
    error('propagateAndCompare:LCPropagationFailed', ...
          'Linked-conic propagation did not reach the 64 SOI boundary.');
end

LC_in = ye_bwd(1, :)';
LC_out = ye_fwd(1, :)';

%% ========================================================================
%  4. N-BODY FULL PROPAGATION AND PERICENTER CONFRONTATION
%  ========================================================================

fprintf('SECTION 4: N-Body Propagation and Pericenter Analysis\n\n');

% --- Setup for the full N-Body propagation from LC_in to LC_out ---
epoch_at_LC_in = pars.INPUTS.epoch0 + te_bwd / 86400;
total_duration_nb = te_fwd - te_bwd; % te_bwd is negative
time_vector_nb = linspace(0, total_duration_nb, pars.GroundTr.npoints)';

% Define the ephemeris function starting from the epoch at LC_in
ephem_func_full_nb = @(t) get_body_positions_wrapper(t, epoch_at_LC_in, actualBodyNaifIDs, spiceParam);

% --- Propagate the full trajectory using the N-Body model ---
[time_out_nb, state_out_nb] = propagateNBodyODE3(LC_in(1:3), LC_in(4:6), ...
                                                 time_vector_nb, ...
                                                 mu_central_body, mu_TBs_local, ...
                                                 ephem_func_full_nb, ...
                                                 specialPerturbationIDs);

% --- Calculate distance to Enceladus for the entire N-body trajectory ---
num_steps_nb = length(time_out_nb);
dist_to_enceladus = zeros(num_steps_nb, 1);
for i = 1:num_steps_nb
    % Get Enceladus position at the corresponding time
    t_sec_from_start_of_nbody = time_out_nb(i);
    r_enc_current = get_body_positions_wrapper(t_sec_from_start_of_nbody, epoch_at_LC_in, pars.INPUTS.NAIFMoon, spiceParam);
    
    % Calculate distance
    dist_to_enceladus(i) = norm(state_out_nb(i, 1:3) - r_enc_current');
end

% --- Find the minimum distance to locate the pericenter ---
[min_dist_nb, idx_pericenter] = min(dist_to_enceladus);

pericenter_NB = state_out_nb(idx_pericenter, :)';
time_at_pericenter_sec_from_NB_start = time_out_nb(idx_pericenter);
epoch_pericenter_NB_mjd = epoch_at_LC_in + time_at_pericenter_sec_from_NB_start / 86400;

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
fprintf('========================================================================\n\n');


%% ========================================================================
%  TEST 1: Forward propagate from LC_in with Kepler to see if we get back
%  ========================================================================

% propagate the same time forward
time_to_propagate_forward = abs(te_bwd);

% Create a new forward time vector for this specific test
time_vector_test_fwd = linspace(eps, time_to_propagate_forward, time_steps)';

% We don't need the SOI crossing event for this test, as we expect to propagate for a known duration.
% So, we'll call the propagator without the event function.
[~, state_test_fwd] = propagateKeplerODE2(LC_in(1:3), LC_in(4:6), time_vector_test_fwd, mu_central_body);

% The final state of this forward propagation is the last row of the state vector
final_state_test = state_test_fwd(end, :)';

% The original initial state
initial_state = [r0_sc_in'; v0_sc_in'];

% --- Compare the results ---
fprintf('TEST 1: Does Kepler forward propagation from LC_in return to the initial state?\n\n');

difference = final_state_test - initial_state;

fprintf('Norm of the difference vector: %e\n', norm(difference));
fprintf('========================================================================\n\n');


%% ========================================================================
%  TEST 2: Forward propagate from LC_in with N-Body to check difference
%  ========================================================================

% The duration is the same as the backward propagation
time_to_propagate_forward_nbody = abs(te_bwd);

% Create the time vector for this test
time_vector_nbody_test = linspace(eps, time_to_propagate_forward_nbody, time_steps)';

% The ephemeris function needs to be aware of the new starting time.
% The backward propagation went from epoch0 to (epoch0 + te_bwd).
% Our new propagation starts at (epoch0 + te_bwd) and goes forward.
epoch_at_LC_in_test = pars.INPUTS.epoch0 + te_bwd / 86400;
ephem_func_nbody_test = @(t) get_body_positions_wrapper(t, epoch_at_LC_in_test, actualBodyNaifIDs, spiceParam);

% --- Propagate forward from LC_in using the N-Body propagator ---
[~, state_nbody_test] = propagateNBodyODE3(LC_in(1:3), LC_in(4:6), ...
                                            time_vector_nbody_test, ...
                                            mu_central_body, mu_TBs_local, ...
                                            ephem_func_nbody_test, ...
                                            specialPerturbationIDs);

% The final state of this N-body forward propagation
final_state_nbody = state_nbody_test(end, :)';

% --- Compare the results with the original initial state ---
fprintf('TEST 2: Forward propagate from LC_in with N-Body propagator.\n\n');

difference_nbody = final_state_nbody - initial_state;

fprintf('Norm of the N-Body difference vector: %e\n', norm(difference_nbody));
fprintf('========================================================================\n\n');


%% ========================================================================
%  7. HELPER FUNCTIONS
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