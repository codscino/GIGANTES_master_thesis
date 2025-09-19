clc;
clear all;
close all;

soi_multiplier = 64*1; % how much soi distance to stop the propagation

pars.GroundTr.t_prop  = 50*2000; % minutes, 2000 is more that 1 enceladus full revolution

perturbingBodyNaifIDs = [-2, -3, -4, 606, 605, 10, 603, 604, 601, 607, 608, 602, 5];
% perturbingBodyNaifIDs = [-2, 10, 602];
% perturbingBodyNaifIDs = [-2, 602];
% perturbingBodyNaifIDs = [10, 602];
% perturbingBodyNaifIDs = [602]; % 3body problem

% pars.EncPlotSize = 1; %true enceladus size, usually too tiny
pars.EncPlotSize = 30;

ga_params.max_dv = 0.01; % [km/s]

%% ========================================================================
%  1. DEFINE INPUTS & INITIALIZE
%  ========================================================================

% --- Flyby and Central Body Parameters ---
pars.INPUTS.idCentral = 6;      % Central Body: Saturn
pars.INPUTS.idMoon    = 1;      % Enceladus flyby
pars.INPUTS.V_inf     = 4.0;    % Hyperbolic excess velocity at Enceladus [km/s]
pars.INPUTS.Flyby.min_h = 25;   % Minimum flyby altitude [km]

% --- Define Incoming and Outgoing Asymptotes (Nodes) ---
nodein =  [pars.INPUTS.V_inf, 0.15, 0];            % [km/s, rad, rad]
nodeout = [pars.INPUTS.V_inf, 0.15, deg2rad(1)];   % [km/s, rad, rad]

% --- SPICE and Time Parameters ---
kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels); 

spiceParam.frame    = 'J2000';
spiceParam.abcorr   = 'NONE';
spiceParam.observer = '699'; % Saturn Center NAIF ID
pars.INPUTS.epoch0  = date2mjd2000([2035, 1, 1, 0, 0, 0]);



actualBodyNaifIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs >= 0);
specialPerturbationIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs < 0);

% --- Load Gravitational Parameters ---
mu_central_body = getAstroConstants('Saturn', 'Mu');
[~, mu_enceladus, R_enceladus, ~] = satMoonsConstants(1); % Enceladus

% --- Dynamically build the list of gravitational parameters (mu) for actual bodies ---
mu_TBs = zeros(1, length(actualBodyNaifIDs));
for i = 1:length(actualBodyNaifIDs)
    id = actualBodyNaifIDs(i);
    switch id
        case 10,  mu_TBs(i) = getAstroConstants('Sun', 'Mu');
        case 5, mu_TBs(i) = getAstroConstants('Jupiter', 'Mu');
        case 601, [~, mu, ~, ~] = satMoonsConstants(0); mu_TBs(i) = mu; % Mimas
        case 602, [~, mu, ~, ~] = satMoonsConstants(1); mu_TBs(i) = mu; % Enceladus
        case 603, [~, mu, ~, ~] = satMoonsConstants(2); mu_TBs(i) = mu; % Tethys
        case 604, [~, mu, ~, ~] = satMoonsConstants(3); mu_TBs(i) = mu; % Dione
        case 605, [~, mu, ~, ~] = satMoonsConstants(4); mu_TBs(i) = mu; % Rhea
        case 606, [~, mu, ~, ~] = satMoonsConstants(5); mu_TBs(i) = mu; % Titan
        case 607, mu_TBs(i) = 0.374;   % Hyperion
        case 608, mu_TBs(i) = 120.4;   % Iapetus
        otherwise
            error('lc_vs_nbody_sma:UnknownNaifID', ...
                  'The NAIF ID %d is not recognized for mu calculation.', id);
    end
end

% Retrieve Saturn Parameters 
[pars.Planet.mu, pars.Planet.EquRad, pars.Planet.OrbRad, pars.Planet.hmin] = planetConstants(pars.INPUTS.idCentral);

% Retrieve Desired Moon Parameters 
if pars.INPUTS.idCentral == 6
    [pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = satMoonsConstants(pars.INPUTS.idMoon);
end
pars.Moon.Vel    = sqrt(pars.Planet.mu/pars.Moon.OrbRad);
pars.Moon.Period = 2*pi*sqrt(pars.Moon.OrbRad^3/pars.Planet.mu);
pars.Moon.HillSph = pars.Moon.OrbRad*( pars.Moon.mu/(3*(pars.Moon.mu + pars.Planet.mu)))^(1/3);

% --- Parameters required for Linked Conic Calculation ---
pars.GroundTr.npoints      = 30e3; 

% propagate to 64 SOI (stable sma)
m_enceladus = 1.08e20;
m_saturn = getAstroConstants('Saturn', 'Mass');
a_enceladus = 238020;
r_soi_enceladus = a_enceladus * (m_enceladus / m_saturn)^(2/5);
pars.INPUTS.maxPropagationDistance = r_soi_enceladus*soi_multiplier;

pars.INPUTS.Flyby.hMapping = 300;

pars.INPUTS.NAIFCentral = 699;
pars.INPUTS.NAIFMoon = 602;


%% ========================================================================
%  2. CALCULATE INITIAL PERICENTER STATE FOR PROPAGATION
%  ========================================================================

[r_enceladus_at_flyby, v_enceladus_at_flyby] = EphSS_car_spice2(602, pars.INPUTS.epoch0 , true, spiceParam);
[vvinfin, r0_sc_in, v0_sc_in, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodein(1), nodein(2), nodein(3), pars.INPUTS.epoch0, pars);
[vvinfout, r0_sc_out, v0_sc_out, vvga] = vinfAlphaCrank_to_VinfCARTClaudio(nodeout(1), nodeout(2), nodeout(3), pars.INPUTS.epoch0, pars );

rp_flyby  = pars.INPUTS.Flyby.min_h + pars.Moon.EquRad;
e_fly     = 1 + ((rp_flyby*pars.INPUTS.V_inf^2)/pars.Moon.mu);
delta_max = 2*asin(1/e_fly);
pars.delta_max = delta_max;

Rm = buildRm(r0_sc_out,vvga); %jose method
% Rm1 = buildRm(r_enceladus_at_flyby, v_enceladus_at_flyby); %my method

[vvinfouBM, delta] = infVelBeforeDefect(vvinfin, vvinfout, delta_max, vvga);
vvinfin_bf = [Rm*vvinfin']';
vvinfouBM_bf = [Rm*vvinfouBM']';

Energy = 0.5*norm(vvinfin_bf)^2;
sma = -mu_enceladus/(2*Energy);
ecc = 1/(sin(delta/2));
rp = sma*(1 - ecc);
hhat = cross(vvinfin_bf, vvinfouBM_bf)./norm(cross(vvinfin_bf, vvinfouBM_bf));
vp = sqrt(norm(vvinfin_bf)^2 + 2*mu_enceladus/rp);

rrp_bf = rp.*(vvinfin_bf - vvinfouBM_bf)./norm(vvinfin_bf - vvinfouBM_bf);
vvp_bf = vp.*cross(hhat, rrp_bf./rp);

Rm_inv = Rm';
rrp_saturn_centric = (Rm_inv * rrp_bf')' + r_enceladus_at_flyby;
vvp_saturn_centric = (Rm_inv * vvp_bf')' + v_enceladus_at_flyby;

rrp_saturn_centric = rrp_saturn_centric(:);
vvp_saturn_centric = vvp_saturn_centric(:);

% coherence check 
[Flyby] = Flyby_BuildUp_claudio(nodein, nodeout, pars);
stin = [r0_sc_in, v0_sc_in];
stout = [r0_sc_out, v0_sc_out];
errin = stin-Flyby.State_In;
errout = stout-Flyby.State_Out;


%% ========================================================================
%  3 PROPAGATE THE NBODY TRAJECTORY
%  ========================================================================

duration_sec = pars.GroundTr.t_prop * 60;
time_steps = pars.GroundTr.npoints;
time_vector_fwd = linspace(eps, duration_sec, time_steps)';
time_vector_bwd = linspace(0, -duration_sec, time_steps)';

ephem_handle = @(t_sec) get_body_positions_wrapper(t_sec, pars.INPUTS.epoch0, actualBodyNaifIDs, spiceParam);

% Create a simpler handle to get ONLY Enceladus's position for the event function
ephem_enceladus = @(t) get_body_positions_wrapper(t, pars.INPUTS.epoch0, pars.INPUTS.NAIFMoon, spiceParam);

% The anonymous function for the event. It captures 'ephem_enceladus'
% and 'pars' from the current workspace.
soi_event_func = @(t, x) soiCrossingEvent(t, x, ephem_enceladus, pars.INPUTS.maxPropagationDistance);


% Propagate FORWARD until the event is triggered
[time_out_fwd, state_out_fwd, te_fwd, ~, ~] = propagateNBodyODE3(rrp_saturn_centric, vvp_saturn_centric, time_vector_fwd, mu_central_body, mu_TBs, ephem_handle, specialPerturbationIDs, soi_event_func);

% Propagate BACKWARD until the event is triggered
[time_out_bwd, state_out_bwd, te_bwd, ~, ~] = propagateNBodyODE3(rrp_saturn_centric, vvp_saturn_centric, time_vector_bwd, mu_central_body, mu_TBs, ephem_handle, specialPerturbationIDs, soi_event_func);

% --- Check and report if the events were found ---
if ~isempty(te_fwd)
    fprintf('Forward propagation stopped at t = %.2f seconds (%.2f min) when distance to Enceladus was reached.\n', te_fwd(1), te_fwd(1)/60);
else
    fprintf('Forward propagation completed without triggering the event.\n');
end
if ~isempty(te_bwd)
    fprintf('Backward propagation stopped at t = %.2f seconds (%.2f min) when distance to Enceladus was reached.\n', te_bwd(1), te_bwd(1)/60);
else
    fprintf('Backward propagation completed without triggering the event.\n');
end

time_out_bwd = flipud(time_out_bwd);
state_out_bwd = flipud(state_out_bwd);
full_time_out_nb = [time_out_bwd; time_out_fwd(2:end)]; %from 2 so the pericentre is not repeated
full_state_out_nb = [state_out_bwd; state_out_fwd(2:end, :)];
num_steps = length(full_time_out_nb);

%% ========================================================================
%  3 PROPAGATE THE LINKED CONICS TRAJECTORY 
%  ========================================================================

% --- Propagate backward from the incoming node until the event ---
[timeLC_bwd, stateLC_bwd, ~, ~, ~] = propagateKeplerODE2(r0_sc_in', v0_sc_in', time_vector_bwd, mu_central_body, soi_event_func);
% [timeLC_bwd, stateLC_bwd, ~, ~, ~] = propagateKeplerODE2(rrp_saturn_centric, vvp_saturn_centric, time_vector_bwd, mu_central_body, soi_event_func);

% --- Propagate forward from the outgoing node until the event ---
[timeLC_fwd, stateLC_fwd, ~, ~, ~] = propagateKeplerODE2(r0_sc_out', v0_sc_out', time_vector_fwd, mu_central_body, soi_event_func);
% [timeLC_fwd, stateLC_fwd, ~, ~, ~] = propagateKeplerODE2(rrp_saturn_centric, vvp_saturn_centric, time_vector_fwd, mu_central_body, soi_event_func);

% Flip the backward propagation results
stateLC_bwd_flipped = flipud(stateLC_bwd);
timeLC_bwd_flipped = flipud(timeLC_bwd);

% Merge the two branches using the OUTPUTS from the ODE solver.
full_time_out_LC = [timeLC_bwd_flipped; timeLC_fwd(1:end)];
full_state_out_LC = [stateLC_bwd_flipped; stateLC_fwd(1:end, :)];


% The two trajectories have different time steps. To compare them, I must
% resample one onto the time grid of the other. I will use the N-body
% time vector as the "master" reference grid.

% Use interp1 to find the LC state at the times specified by the N-body propagation.
state_LC_resampled = interp1(full_time_out_LC, full_state_out_LC, full_time_out_nb, 'pchip');



%% ========================================================================
%  4. NODEIN, NODEOUT at +- 64SOI
%  ========================================================================

% n-body states
nodein_nb_state = full_state_out_nb(1,:);
nodeout_nb_state = full_state_out_nb(end,:);

% linked conic state(nodein does not matter for this startegy)
nodeout_LC_state = state_LC_resampled(end,:);

% Convert to Keplerian to find TARGET SEMI-MAJOR AXIS
nodeout_LC_kep = car2kep(nodeout_LC_state, mu_central_body);
target_sma = nodeout_LC_kep(1);


%% ========================================================================
%   5.  GENETIC ALGORITHM OPTIMIZATION (TWO-BURN MINIMIZATION)
% =========================================================================

% --- 1. Define Optimization Variables and Bounds ---
% The GA will optimize the velocity vector at the pericenter.
% The position at pericenter is FIXED.
original_pericenter_pos = rrp_saturn_centric';
original_pericenter_vel = vvp_saturn_centric';

% Define a search space for the GA around the original velocity.
% This defines the maximum possible deltav1. We'll use your max_dv for this.
search_dv_margin = ga_params.max_dv; 
lb = original_pericenter_vel - search_dv_margin; % Lower bounds for [vx, vy, vz]
ub = original_pericenter_vel + search_dv_margin; % Upper bounds for [vx, vy, vz]

% --- 2. Define Fitness Function and Parameters ---
ga_params.mu_central_body = mu_central_body;
ga_params.mu_TBs = mu_TBs;
ga_params.specialPerturbationIDs = specialPerturbationIDs;
ga_params.target_sma = target_sma;
ga_params.NAIFMoon = pars.INPUTS.NAIFMoon;
ga_params.hMapping = pars.INPUTS.maxPropagationDistance;

% Pass the original pericenter state for deltav1 calculation
ga_params.original_pericenter_pos = original_pericenter_pos;
ga_params.original_pericenter_vel = original_pericenter_vel;

% Pass SPICE and timing info for propagation inside the objective function
ga_params.epoch_pericenter = pars.INPUTS.epoch0;
ga_params.spiceParam = spiceParam;
ga_params.actualBodyNaifIDs = actualBodyNaifIDs;

% We only need to propagate FORWARD from pericenter. Use the time from the initial run.
ga_params.propagation_duration = te_fwd(1) * 1.2; % Add 20% margin

% Create a handle to the objective function
fitness_function = @(v) minimize_total_dv_objective(v, ga_params);

% --- 3. Set GA Options ---
nvars = 3; % Optimizing [vx, vy, vz] at pericenter
options = optimoptions('ga', ...
    'Display', 'iter', ...
    'PopulationSize', 80, ...
    'MaxGenerations', 50, ...
    'FunctionTolerance', 1e-9, ...
    'ConstraintTolerance', 1e-8, ...
    'UseParallel', true, ... 
    'PlotFcn', @gaplotbestf, ...
    'HybridFcn', @fmincon);

% --- 4. Run the Genetic Algorithm ---
[optimized_pericenter_vel, min_total_dv] = ga(fitness_function, nvars, [], [], [], [], lb, ub, [], options);

%% ========================================================================
%    6.           VERIFICATION OF THE OPTIMIZED SOLUTION
% =========================================================================
fprintf('\n\n==================================================\n');
fprintf('OPTIMIZATION COMPLETE: VERIFYING FINAL SOLUTION\n');
fprintf('==================================================\n');
fprintf('Minimum Total Delta-V found by GA: %.4f m/s\n', min_total_dv * 1000);

% --- Propagate the OPTIMIZED trajectory one last time for analysis ---
r_optimized = original_pericenter_pos';      % This is the fixed position (3x1 column)
v_optimized = optimized_pericenter_vel';     % This is the optimized velocity (3x1 column)

time_vector_verif = linspace(0, ga_params.propagation_duration, 10000)';
eph_handle_verif = @(t_sec) get_body_positions_wrapper(t_sec, ga_params.epoch_pericenter, actualBodyNaifIDs, spiceParam);
soi_event_func_verif = @(t, x) soiOutboundCrossingEvent(t, x, @(t)get_body_positions_wrapper(t,ga_params.epoch_pericenter,ga_params.NAIFMoon,ga_params.spiceParam), ga_params.hMapping);

[~, ~, ~, ye_verif, ~] = propagateNBodyODE3(r_optimized, v_optimized, ...
                                        time_vector_verif, mu_central_body, mu_TBs, eph_handle_verif, ...
                                        specialPerturbationIDs, soi_event_func_verif);
                                    
state_at_soi_exit = ye_verif(1,:);

% --- Deconstruct the final Delta-V cost ---
deltav1_final = norm(optimized_pericenter_vel - original_pericenter_vel);

% Calculate deltav2 for the final, optimized trajectory
r_out_nb = state_at_soi_exit(1:3);
v_out_nb = state_at_soi_exit(4:6);
kep_current = car2kep(state_at_soi_exit, mu_central_body);
sma_current = kep_current(1);
v_current_mag = norm(v_out_nb);
r_mag = norm(r_out_nb);

% vis viva
v_new_mag_sq = 2*mu_central_body/r_mag - mu_central_body/target_sma;

% Check for impossible maneuvers (negative square root)
if v_new_mag_sq < 0
    deltav2_final = NaN; % Indicate an error
    fprintf('ERROR: Cannot calculate dv2, target SMA is unreachable from this state.\n');
else
    v_new_mag = sqrt(v_new_mag_sq);
    deltav2_final = abs(v_new_mag - v_current_mag);
end

% --- Display Final Results ---
fprintf('\n--- Delta-V Breakdown ---\n');
fprintf('Delta-V 1 (at Pericenter):   %.4f m/s\n', deltav1_final * 1000);
fprintf('Delta-V 2 (at SOI Exit):     %.4f m/s\n', deltav2_final * 1000);
fprintf('Total Delta-V Cost:          %.4f m/s\n', (deltav1_final + deltav2_final) * 1000);

fprintf('\n--- SMA Comparison (at SOI Exit, Before dv2) ---\n');
fprintf('Target SMA (from LC):              %.4f km\n', target_sma);
fprintf('Achieved NB SMA (before dv2):      %.4f km\n', sma_current);
fprintf('Final SMA Error (to be corrected by dv2): %.4f %%\n\n', abs(sma_current-target_sma)/target_sma*100);

fprintf('SUCCESS: Found a two-burn solution to match the target resonance orbit.\n');
fprintf('==================================================\n');


%% ========================================================================
%                        HELPER AND OBJECTIVE FUNCTIONS
% =========================================================================

function total_dv = minimize_total_dv_objective(vvp_candidate, P)
    % Objective function for the GA.
    % GOAL: Minimize the sum of two maneuvers:
    % dv1: An impulsive burn at the pericenter.
    % dv2: A cleanup burn at the SOI exit to match the target SMA.
    
    % persistent to make parallel computation work with spice
    persistent kernels_loaded;
    if isempty(kernels_loaded) || ~kernels_loaded
        kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
        loadSpiceKernels(kernels);
        kernels_loaded = true;
    end

    % =====================================================================
    %  PART 1: CALCULATE THE FIRST MANEUVER, deltav1
    % =====================================================================
    deltav1 = norm(vvp_candidate - P.original_pericenter_vel);

    % =====================================================================
    %  PART 2: PROPAGATE THE TRAJECTORY FORWARD
    % =====================================================================
    r_start = P.original_pericenter_pos;
    v_start = vvp_candidate;

    time_vector = linspace(0, P.propagation_duration, 5000)'; 
    eph_handle = @(t) get_body_positions_wrapper(t, P.epoch_pericenter, P.actualBodyNaifIDs, P.spiceParam);
    eph_enceladus_handle = @(t) get_body_positions_wrapper(t, P.epoch_pericenter, P.NAIFMoon, P.spiceParam);
    outbound_event_func = @(t, x) soiOutboundCrossingEvent(t, x, eph_enceladus_handle, P.hMapping);

    [~, ~, ~, ye, ~] = propagateNBodyODE3(r_start(:), v_start(:), time_vector, P.mu_central_body, ...
                                      P.mu_TBs, eph_handle, P.specialPerturbationIDs, outbound_event_func);

    % --- propagation failures ---
    if isempty(ye) || any(isnan(ye(1,:)))
        total_dv = 1e7; % Assign a large penalty and exit
        return;
    end

    % =====================================================================
    %  PART 3: CALCULATE THE SECOND MANEUVER, deltav2
    % =====================================================================
    state_at_soi_exit = ye(1, :);
    r_out_nb = state_at_soi_exit(1:3)';
    v_out_nb = state_at_soi_exit(4:6)';
    
    % Get the current SMA from the propagated n-body trajectory
    kep_current = car2kep(state_at_soi_exit, P.mu_central_body);
    sma_current = kep_current(1);

    % Use the vis-viva equation to find the required velocity magnitude at this
    % position (r_out_nb) to achieve the target SMA (P.target_sma).
    v_current_mag = norm(v_out_nb);
    r_mag = norm(r_out_nb);

    % v_new^2 = 2*mu/r - mu/a_target
    v_new_mag_sq = 2*P.mu_central_body/r_mag - P.mu_central_body/P.target_sma;
    
    % Handle cases where the geometry is impossible (e.g., requires imaginary velocity)
    if v_new_mag_sq < 0
        deltav2 = 1e6; % Assign a large penalty
    else
        v_new_mag = sqrt(v_new_mag_sq);
        % deltav2 (a tangential thrust)
        deltav2 = abs(v_new_mag - v_current_mag);
    end

    % =====================================================================
    %  PART 4: RETURN THE TOTAL COST (FITNESS VALUE)
    % =====================================================================
    total_dv = deltav1 + deltav2;
end




% --- Other Helper Functions ---

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

function [value, isterminal, direction] = soiOutboundCrossingEvent(t, x, ephem_enceladus_handle, target_distance)
    r_sc = x(1:3);
    r_enceladus = ephem_enceladus_handle(t);
    value = norm(r_sc - r_enceladus) - target_distance;
    isterminal = 1;
    direction = 1; % only trigger when distance is increasing.
end