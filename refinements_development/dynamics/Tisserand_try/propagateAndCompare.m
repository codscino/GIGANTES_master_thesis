function [pericenter_LC, pericenter_NB, epoch_pericenter_NB, LC_out, NB_out, T_in, T_out] = ...
    propagateAndCompare(nodein, nodeout, pars)

% propagates and compares linked-conic and n-body trajectories.
%
% This function takes nodein and nodeout Enceladus flyby paramters and all
% the costants + the epoch of the time of the pericentre from pars
% structure
%
% to perform the following steps:
% 1.  Calculates the pericenter state of a flyby trajectory using linked-conics.
% 2.  Propagates from Enceladus center LC 2 body out to +/- 64 SOI of Enceladus (LC_in, LC_out).
% 3.  Propagates in n-body from NB_in = LC_in forward to +64 SOI (NB_out).
% 4.  Determines the pericenter of the n-body trajectory.
% 5.  Calculates the Tisserand parameter at the start and end of the n-body propagation.
%
% OUTPUTS:
%   pericenter_LC     - Cartesian state [pos; vel] of the LC pericenter (saturn centric column vector)
%   pericenter_NB     - Cartesian state [pos; vel] of the NB pericenter. (saturn centric column vector)
%   epoch_pericenter_NB - Epoch (MJD) of the n-body pericenter (days)
%   LC_out            - Cartesian state of the linked-conic trajectory at +64 SOI. (saturn centric column vector)
%   NB_out            - Cartesian state of the n-body trajectory at +64 SOI. (saturn centric column vector)
%   T_in              - Tisserand parameter with respect to Enceladus at the start (NB_in).
%   T_out             - Tisserand parameter with respect to Enceladus at the end (NB_out).


%% ========================================================================
%  1. INITIALIZE PARAMETERS & SPICE
%  ========================================================================

spiceParam.frame    = 'J2000';
spiceParam.abcorr   = 'NONE';
spiceParam.observer = '699'; % Saturn Center NAIF ID

mu_central_body = pars.Planet.mu;
mu_enceladus = pars.Moon.mu;

% Perturbing bodies setup from the 'pars' struct
perturbingBodyNaifIDs = pars.INPUTS.perturbingBodyNaifIDs;
actualBodyNaifIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs >= 0);
specialPerturbationIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs < 0);
mu_TBs = pars.INPUTS.mu_TBs;


%% ========================================================================
%  2. CALCULATE LINKED-CONIC PERICENTER AND IN/OUT STATES
%  ========================================================================

% --- Get Enceladus' state at the flyby epoch ---
[r_enceladus_at_flyby, v_enceladus_at_flyby] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, pars.INPUTS.epoch0, true, spiceParam);

% --- Calculate incoming and outgoing V-infinity vectors ---
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

pericenter_LC = [rrp_saturn_centric(:); vvp_saturn_centric(:)];

%% ========================================================================
%  3. PROPAGATE LINKED CONICS TO +/- 64 ENCELADUS SOI
%  ========================================================================

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
[timeLC_bwd, stateLC_bwd, te_bwd, ye_bwd, ~] = propagateKeplerODE2(r0_sc_in', v0_sc_in', time_vector_bwd, mu_central_body, soi_event_func);
[timeLC_fwd, stateLC_fwd, te_fwd, ye_fwd, ~] = propagateKeplerODE2(r0_sc_out', v0_sc_out', time_vector_fwd, mu_central_body, soi_event_func);

if isempty(ye_fwd) || isempty(ye_bwd)
    error('propagateAndCompare:LCPropagationFailed', ...
          'Linked-conic propagation did not reach the 64 SOI boundary.');
end

LC_in = ye_bwd(1, :)';
LC_out = ye_fwd(1, :)';

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

% Use the new event function to detect the EXIT boundary crossing only
soi_event_func_nbody = @(t, x) soiExitEvent(t, x, ephem_enceladus_nbody, soi_64_dist); % <-- CHANGE HERE

[time_out_nb, state_out_nb, te_fwd_nb, ye_fwd_nb, ~] = propagateNBodyODE3(NB_in(1:3), NB_in(4:6), ...
    time_vector_nbody, mu_central_body, mu_TBs, ephem_handle, specialPerturbationIDs, soi_event_func_nbody);


if isempty(te_fwd_nb)
   warning('propagateAndCompare:NBPropagationWarning', ...
           'N-body propagation completed the full duration without reaching the 64 SOI boundary.');
   NB_out = state_out_nb(end,:)';
else
   NB_out = ye_fwd_nb(1, :)';
end


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

pericenter_NB = state_out_nb(idx_pericenter, :)';
time_at_pericenter_sec_from_NB_start = time_out_nb(idx_pericenter);
epoch_pericenter_NB = initial_epoch_offset_days + time_at_pericenter_sec_from_NB_start / 86400;


%% ========================================================================
%  6. CALCULATE TISSERAND PARAMETERS
%  ========================================================================
kep_in = car2kep(NB_in', mu_central_body);
kep_out = car2kep(NB_out', mu_central_body);

% --- Get Enceladus state at the NB_in epoch ---
[r_enc_in, v_enc_in] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, initial_epoch_offset_days, true, spiceParam);
state_enc_in = [r_enc_in, v_enc_in];
kep_enc_in = car2kep(state_enc_in, mu_central_body);
a_enc_in = kep_enc_in(1);

% --- Get Enceladus state at the NB_out epoch ---
epoch_out = initial_epoch_offset_days + te_fwd_nb / 86400;
[r_enc_out, v_enc_out] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, epoch_out, true, spiceParam);
state_enc_out = [r_enc_out, v_enc_out];
kep_enc_out = car2kep(state_enc_out, mu_central_body);
a_enc_out = kep_enc_out(1);

% --- Calculate Tisserand Parameter at in and out points ---
T_in = a_enc_in/kep_in(1) + 2*cos(kep_in(3)) * sqrt( (kep_in(1)/a_enc_in) * (1- kep_in(2)^2) );
T_out = a_enc_out/kep_out(1) + 2*cos(kep_out(3)) * sqrt( (kep_out(1)/a_enc_out) * (1- kep_out(2)^2) );


end % end of propagateAndCompare


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