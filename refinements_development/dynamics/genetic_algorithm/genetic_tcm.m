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

ga_params.max_dv = 0.003; % [km/s]

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
%  3.1 PROPAGATE THE NBODY TRAJECTORY
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
%  3.2 PROPAGATE THE LINKED CONICS TRAJECTORY 
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

%% ========================================================================
%  3.3 RESAMPLE TRAJECTORIES TO A COMMON TIME GRID FOR COMPARISON
%  ========================================================================

% The two trajectories have different time steps. To compare them, I must
% resample one onto the time grid of the other. I will use the N-body
% time vector as the "master" reference grid.

% Use interp1 to find the LC state at the times specified by the N-body propagation.
state_LC_resampled = interp1(full_time_out_LC, full_state_out_LC, full_time_out_nb, 'pchip');


%% ========================================================================
%  8. NODEIN, NODEOUT at +- 64SOI
%  ========================================================================
nodein_nb = car2kep(full_state_out_nb(1,:), mu_central_body);
nodein_LC = car2kep(state_LC_resampled(1,:), mu_central_body);

nodeout_nb = car2kep(full_state_out_nb(end,:), mu_central_body);
nodeout_LC = car2kep(state_LC_resampled(end,:), mu_central_body);


%% ========================================================================
%       GENETIC ALGORITHM OPTIMIZATION (VELOCITY VECTOR)
% =========================================================================
fprintf('\n\n==================================================\n');
fprintf('STARTING GENETIC ALGORITHM (VELOCITY OPTIMIZATION)\n');
fprintf('==================================================\n');

% --- 1. Define Target, Initial State, and Error ---
target_sma = nodeout_LC(1);
nodeout_nb_initial = car2kep(full_state_out_nb(end,:), mu_central_body);
initial_sma_nb = nodeout_nb_initial(1);
initial_error_percent = abs((initial_sma_nb - target_sma) / target_sma) * 100;

fprintf('Target SMA (from Linked Conics):      %.4f km\n', target_sma);
fprintf('Original N-Body SMA (from propagation): %.4f km\n', initial_sma_nb);
fprintf('Initial SMA difference:               %.4f %%\n\n', initial_error_percent);

% --- 2. Define Optimization Variables and Bounds ---
% The maneuver is applied at the start of the backward propagation.
% The position is FIXED, we only optimize the velocity vector.
original_state_cart = full_state_out_nb(1,:);
original_position_cart = original_state_cart(1:3);
original_velocity_cart = original_state_cart(4:6);

% Define a search space around the original velocity vector.
% The search space should be larger than the max_dv constraint.
search_dv_margin = 1.5 * ga_params.max_dv; 
lb = original_velocity_cart - search_dv_margin; % Lower bounds for [vx, vy, vz]
ub = original_velocity_cart + search_dv_margin; % Upper bounds for [vx, vy, vz]

% --- 3. Define Fitness Function and Parameters ---
ga_params.mu_central_body = mu_central_body;
ga_params.mu_TBs = mu_TBs;
ga_params.specialPerturbationIDs = specialPerturbationIDs;
ga_params.target_sma = target_sma;
ga_params.NAIFMoon = pars.INPUTS.NAIFMoon;
ga_params.hMapping = pars.INPUTS.maxPropagationDistance;

% Pass the FIXED position and ORIGINAL velocity for dV calculation
ga_params.original_position_cart = original_position_cart;
ga_params.original_velocity_cart = original_velocity_cart;

% The propagation starts at the time of the backward SOI crossing.
ga_params.epoch_start_prop = pars.INPUTS.epoch0 + te_bwd(1) / 86400;
ga_params.spiceParam = spiceParam;
ga_params.actualBodyNaifIDs = actualBodyNaifIDs;

% Total duration of the flyby propagation. A margin is added.
ga_params.propagation_duration = (te_fwd(1) - te_bwd(1)) * 1.2;

% Create a handle to the NEW objective function
fitness_function = @(v) sma_error_objective_velocity(v, ga_params);

% --- 4. Set GA Options ---
nvars = 3; % We are only optimizing [vx, vy, vz]
options = optimoptions('ga', ...
    'Display', 'iter', ...
    'PopulationSize', 60, ...
    'MaxGenerations', 150, ...
    'FunctionTolerance', 1e-8, ...
    'ConstraintTolerance', 1e-7, ...
    'UseParallel', true, ... 
    'PlotFcn', @gaplotbestf);

% --- 5. Run the Genetic Algorithm ---
fprintf('Running Genetic Algorithm... Optimizing initial velocity vector.\n\n');
[optimized_velocity, fval] = ga(fitness_function, nvars, [], [], [], [], lb, ub, [], options);

%% ========================================================================
%               VERIFICATION OF THE OPTIMIZED SOLUTION
% =========================================================================
fprintf('\n\n==================================================\n');
fprintf('OPTIMIZATION COMPLETE: VERIFYING FINAL SOLUTION\n');
fprintf('==================================================\n');
fprintf('Final objective function value (percentage error): %.8f\n', fval*100);

% --- Propagate with the optimized initial state for verification ---
% The initial state for propagation uses the FIXED original position
% and the NEWLY optimized velocity vector.
state_in_optimized_cart = [original_position_cart, optimized_velocity];

% Define time vector for the final, full-resolution verification
time_vector_verif = linspace(0, ga_params.propagation_duration, time_steps)';
eph_handle_verif = @(t_sec) get_body_positions_wrapper(t_sec, ga_params.epoch_start_prop, actualBodyNaifIDs, spiceParam);

% Propagate the corrected trajectory
[~, state_out_verif] = propagateNBodyODE3(state_in_optimized_cart(1:3)', state_in_optimized_cart(4:6)', ...
                                        time_vector_verif, mu_central_body, mu_TBs, eph_handle_verif, ...
                                        specialPerturbationIDs, []); % No event function needed for verification

% Get the final state and convert to Keplerian elements for comparison
final_state_cart_optimized = state_out_verif(end, :);
final_node_kep_optimized = car2kep(final_state_cart_optimized, mu_central_body);

% --- Compare results ---
final_sma_optimized = final_node_kep_optimized(1);
final_error_percent = abs((final_sma_optimized - target_sma) / target_sma) * 100;

fprintf('\n--- Velocity Comparison (at maneuver point) ---\n');
fprintf('Original Velocity:   [%.4f, %.4f, %.4f] km/s\n', original_velocity_cart);
fprintf('Optimized Velocity:  [%.4f, %.4f, %.4f] km/s\n', optimized_velocity);

fprintf('\n--- SMA Comparison ---\n');
fprintf('Target SMA (from LC):            %.4f km\n', target_sma);
fprintf('Final NB SMA (from optimized):   %.4f km\n', final_sma_optimized);
fprintf('Final Percentage Difference:     %.6f %%\n\n', final_error_percent);

% --- Calculate and Display the Delta-V for the Final Solution ---
% The Delta-V is now a simple, direct vector norm calculation.
fprintf('--- Delta-V Cost for Optimized Solution ---\n');
final_dv = norm(optimized_velocity - original_velocity_cart);

fprintf('Minimum Delta-V to correct trajectory: %.4f m/s\n', final_dv * 1000);
fprintf('Constraint was max_dv = %.4f m/s\n', ga_params.max_dv * 1000);

if final_error_percent <= 0.1 && final_dv <= ga_params.max_dv
    fprintf('\nSUCCESS: The optimized semi-major axis is within the target tolerance AND the dV is within budget.\n');
else
    fprintf('\nFAILURE: The optimization did not meet all constraints.\n');
end
fprintf('==================================================\n');

%% ========================================================================
%                        HELPER AND OBJECTIVE FUNCTIONS
% =========================================================================

function error = sma_error_objective_velocity(v_candidate, P)
    % Objective function for the genetic algorithm (Velocity Optimization).
    % 1. Penalizes solutions with delta-V > max_dv.
    % 2. Propagates from a fixed position with the candidate velocity.
    % 3. Minimizes the error in the final semi-major axis.
    persistent kernels_loaded;
    if isempty(kernels_loaded) || ~kernels_loaded
        kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
        loadSpiceKernels(kernels);
        kernels_loaded = true;
    end

    % =====================================================================
    %  PART 1: DELTA-V CONSTRAINT CHECK (IMPULSIVE MANEUVER)
    % =====================================================================
    % Calculate the dV cost as a single impulsive burn
    delta_v = norm(v_candidate - P.original_velocity_cart);

    % --- Apply Penalty if dV is too high ---
    if delta_v > P.max_dv
        % Assign a large penalty. The penalty increases the further it is
        % from the valid dV, which can help guide the search.
        error = 1e6 + (delta_v - P.max_dv)^2;
        return;
    end

    % =====================================================================
    %  PART 2: N-BODY PROPAGATION WITH OUTBOUND EVENT
    % =====================================================================
    % The starting position is fixed. The GA only provides the velocity.
    r_start = P.original_position_cart';
    v_start = v_candidate';

    % --- Define propagation settings ---
    time_vector = linspace(0, P.propagation_duration, 2000)';
    eph_handle = @(t) get_body_positions_wrapper(t, P.epoch_start_prop, P.actualBodyNaifIDs, P.spiceParam);
    eph_enceladus = @(t) get_body_positions_wrapper(t, P.epoch_start_prop, P.NAIFMoon, P.spiceParam);
    outbound_event_func = @(t, x) soiOutboundCrossingEvent(t, x, eph_enceladus, P.hMapping);

    % --- Propagate the trajectory until the outbound event triggers ---
    [~, ~, te, ye, ~] = propagateNBodyODE3(r_start, v_start, time_vector, P.mu_central_body, ...
                                      P.mu_TBs, eph_handle, P.specialPerturbationIDs, outbound_event_func);

    % --- Handle propagation failures or if the event was not found ---
    if isempty(te) || any(isnan(ye(1,:)))
        error = 5e5; % Assign a large penalty if the flyby fails to complete
        return;
    end

    % --- Get the final state at the event and convert to Keplerian ---
    final_state_cart = ye(1, :);
    final_kep = car2kep(final_state_cart, P.mu_central_body);
    final_sma = final_kep(1);

    % --- Calculate the final error (fitness value) ---
    error = abs((final_sma - P.target_sma) / P.target_sma);

    if isnan(error)
        error = 1e6; % Assign a penalty for any other NaN results
    end
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
    direction = 1; % Key difference: only trigger when distance is increasing.
end