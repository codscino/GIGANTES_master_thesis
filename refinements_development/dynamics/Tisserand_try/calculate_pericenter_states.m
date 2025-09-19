function [pLC, pNB] = calculate_pericenter_states(nodein, nodeout, pars)
% calculate_pericenter_states: Computes and compares pericenter states.
%
% This function calculates the pericenter state of a spacecraft during a flyby
% of a moon (e.g., Enceladus) using two different models: a linked-conics
% approximation and a full n-body numerical propagation.
%
% INPUTS:
%   nodein  (1x3 vector): Incoming V-infinity vector components 
%                         [V_inf, B-plane angle, crank angle].
%   nodeout (1x3 vector): Outgoing V-infinity vector components.
%   pars    (struct):     A structure containing all necessary simulation
%                         parameters, including planetary constants, SPICE
%                         kernel info, and propagation settings.
%
% OUTPUTS:
%   pLC (6x1 vector): The state vector [r; v] at the pericenter as
%                     calculated by the linked-conics model.
%   pNB (6x1 vector): The state vector [r; v] at the pericenter as found
%                     from the n-body propagation.
%

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
% NOTE: Replaced 'nodein_base' and 'nodeout_base' with function inputs
[vvinfin, r0_sc_in, v0_sc_in, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodein(1), nodein(2), nodein(3), pars.INPUTS.epoch0, pars);
[vvinfout, r0_sc_out, v0_sc_out, vvga] = vinfAlphaCrank_to_VinfCARTClaudio(nodeout(1), nodeout(2), nodeout(3), pars.INPUTS.epoch0, pars);

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

% Assign to the first output variable
pLC = [rrp_saturn_centric(:); vvp_saturn_centric(:)];


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

%% Section 4: N-BODY FULL PROPAGATION AND PERICENTER CONFRONTATION
%  --------------------------------------------------------------------

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
[~, idx_pericenter] = min(dist_to_enceladus);

% Assign to the second output variable
pNB = state_out_nb(idx_pericenter, :)';

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
        positions(:, k) = r(:); % Ensure column vector
    end
end

function [value, isterminal, direction] = soiCrossingEvent(t, x, ephem_enceladus_handle, target_distance)
    r_sc = x(1:3);
    r_enceladus = ephem_enceladus_handle(t);
    current_distance = norm(r_sc - r_enceladus);
    value = current_distance - target_distance;
    isterminal = 1; % Stop the integration
    direction = 0;  % Detect crossing from any direction
end