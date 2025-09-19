%% STM with also deltaV calculation, and optimization
% the otimization is multi varibale(B, deltaV) using least squares and
% variable weights

clc;
clear all;
close all;

%% ========================================================================
%  1. SETUP 
%  ========================================================================
soi_multiplier = 640; 

% Your existing parameter setup
pars.t_prop_hours = 240;
pars.GroundTr.npoints = 30e3;
pars.INPUTS.perturbingBodyNaifIDs = [];
% pars.INPUTS.perturbingBodyNaifIDs = [602,10];
% pars.INPUTS.perturbingBodyNaifIDs = [-2, 10];
% pars.INPUTS.perturbingBodyNaifIDs = [-2, 602, 10]; % J2, Enceladus, Sun
pars.INPUTS.idCentral = 6;
pars.INPUTS.idMoon = 1;
pars.INPUTS.Flyby.min_h = 10;

% Load kernels and constants (your existing code)
kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels);

[pars.Planet.mu, ~, ~, ~] = planetConstants(pars.INPUTS.idCentral);
[pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = satMoonsConstants(pars.INPUTS.idMoon);
pars.Moon.HillSph = pars.Moon.OrbRad * (pars.Moon.mu / (3 * pars.Planet.mu))^(1/3);

% Your other setup parameters
pars.INPUTS.NAIFCentral = 699;
pars.INPUTS.NAIFMoon = 602;


% propagate to 64 SOI (stable sma 3 body)
m_enceladus = 1.08e20;
m_saturn = getAstroConstants('Saturn', 'Mass');
a_enceladus = 238020;
r_soi_enceladus = a_enceladus * (m_enceladus / m_saturn)^(2/5);

pars.INPUTS.maxPropagationDistance = r_soi_enceladus*soi_multiplier;

pars.INPUTS.epoch0 = date2mjd2000([2040 1 1 12 0 0]);
pars.INPUTS.V_inf = 4;

% Flyby parameters
nodein = [pars.INPUTS.V_inf, 0, 0];
nodeout = [pars.INPUTS.V_inf, 0.001, deg2rad(1)];


%% ========================================================================
%  3. RUN MULTI-OBJECTIVE STM OPTIMIZATION
%  ========================================================================

fprintf('\n========================================\n');
fprintf('MULTI-OBJECTIVE STM OPTIMIZATION\n');
fprintf('========================================\n');

% Call the multi-objective optimization function
[NB_in_corrected, B_achieved, deltaV_achieved] = performSTMOptimizationMultiObjective(pars, nodein, nodeout);

% Check if objectives were met
if deltaV_achieved < 20  % 20 m/s tolerance
    fprintf('\n✓ SUCCESS: Delta-V within tolerance (%.3f m/s)\n', deltaV_achieved);
else
    fprintf('\n⚠ WARNING: Delta-V exceeds tolerance (%.3f m/s > 20 m/s)\n', deltaV_achieved);
    fprintf('Consider adjusting nodeout parameters or accepting higher delta-V.\n');
end

%% ========================================================================
%  4. VERIFY THE CORRECTED TRAJECTORY
%  ========================================================================

fprintf('\nVerifying corrected trajectory...\n');

% Setup for verification propagation
soi_multiplier = 64;
spiceParam.frame = 'J2000';
spiceParam.abcorr = 'NONE';
spiceParam.observer = '699';

% Get initial epoch at -64 SOI
duration_sec = pars.t_prop_hours * 3600;
time_steps = pars.GroundTr.npoints;

% Calculate when we're at -64 SOI (you already have this from your backpropagation)
[r_enceladus_at_flyby, v_enceladus_at_flyby] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, pars.INPUTS.epoch0, true, spiceParam);
[vvinfin, r0_sc_in, v0_sc_in, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodein(1), nodein(2), nodein(3), pars.INPUTS.epoch0, pars);

time_vector_bwd = linspace(0, -duration_sec, time_steps)';
soi_64_dist = pars.INPUTS.maxPropagationDistance;
ephem_enceladus = @(t) get_body_positions_wrapper(t, pars.INPUTS.epoch0, pars.INPUTS.NAIFMoon, spiceParam);
soi_event_func = @(t, x) soiCrossingEvent(t, x, ephem_enceladus, soi_64_dist);

[~, ~, te_bwd, ~, ~] = propagateKeplerODE2(r0_sc_in', v0_sc_in', time_vector_bwd, pars.Planet.mu, soi_event_func);
initial_epoch_offset_days = pars.INPUTS.epoch0 + te_bwd/86400;

% Propagate the corrected trajectory
flyby_tot_duration = duration_sec * 2;
time_vector_nbody = linspace(eps, flyby_tot_duration, time_steps)';

% Build perturber lists
perturbingBodyNaifIDs = pars.INPUTS.perturbingBodyNaifIDs;
specialPerturbationIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs < 0);
actualBodyNaifIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs >= 0);

% Build mu_TBs
mu_TBs = [];
for i = 1:length(actualBodyNaifIDs)
    id = actualBodyNaifIDs(i);
    switch id
        case 10,  mu_TBs(i) = getAstroConstants('Sun', 'Mu');
        case 602, [~, mu, ~, ~] = satMoonsConstants(1); mu_TBs(i) = mu;
        % Add other cases as needed
    end
end

ephem_handle = @(t_sec) get_body_positions_wrapper(t_sec, initial_epoch_offset_days, actualBodyNaifIDs, spiceParam);
ephem_enceladus_nbody = @(t) get_body_positions_wrapper(t, initial_epoch_offset_days, pars.INPUTS.NAIFMoon, spiceParam);

[time_out_verify, state_out_verify, ~, ~, ~] = propagateNBodyODE3(NB_in_corrected(1:3), NB_in_corrected(4:6), ...
    time_vector_nbody, pars.Planet.mu, mu_TBs, ephem_handle, specialPerturbationIDs, []);

% Find periapsis in verification
num_steps = length(time_out_verify);
dist_to_enceladus = zeros(num_steps, 1);
for i = 1:num_steps
    r_enc_current = ephem_enceladus_nbody(time_out_verify(i));
    dist_to_enceladus(i) = norm(state_out_verify(i, 1:3) - r_enc_current');
end

[min_dist_verify, idx_peri_verify] = min(dist_to_enceladus);
altitude_verify = min_dist_verify - pars.Moon.EquRad;

fprintf('\nVerification Results:\n');
fprintf('--------------------\n');
fprintf('Achieved periapsis altitude: %.3f km\n', altitude_verify);
fprintf('Minimum distance to Enceladus: %.3f km\n', min_dist_verify);

% Calculate final B-plane
state_at_peri_verify = state_out_verify(idx_peri_verify, :)';
time_at_peri_verify = time_out_verify(idx_peri_verify);
epoch_peri_verify = initial_epoch_offset_days + time_at_peri_verify/86400;

[r_enc_peri, v_enc_peri] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, epoch_peri_verify, true, spiceParam);
r_sc_rel = state_at_peri_verify(1:3) - r_enc_peri';
v_sc_rel = state_at_peri_verify(4:6) - v_enc_peri';

v_inf_sq = dot(v_sc_rel, v_sc_rel) - 2*pars.Moon.mu/norm(r_sc_rel);
[B_R_verify, B_T_verify, ~] = b_plane_targeting(r_sc_rel, v_sc_rel, v_inf_sq, pars.Moon.mu);

fprintf('\nFinal B-plane parameters:\n');
fprintf('  B_R: %.6f km\n', B_R_verify);
fprintf('  B_T: %.6f km\n', B_T_verify);
fprintf('  |B|: %.6f km\n', norm([B_R_verify; B_T_verify]));

%% ========================================================================
%  5. COMPARISON PLOT
%  ========================================================================

% add visualization here to compare:
% - Original linked conics trajectory
% - Uncorrected N-body trajectory  
% - STM-corrected N-body trajectory

fprintf('\n========================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('========================================\n');

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

function [value, isterminal, direction] = soiCrossingEvent(t, x, ephem_enceladus_handle, target_distance)
    r_sc = x(1:3);
    r_enceladus = ephem_enceladus_handle(t);
    current_distance = norm(r_sc - r_enceladus);
    value = current_distance - target_distance;
    isterminal = 1;
    direction = 0;
end