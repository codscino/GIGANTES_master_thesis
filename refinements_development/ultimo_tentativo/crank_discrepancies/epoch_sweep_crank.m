clc;
clear all;
close all;

pars.GroundTr.t_prop  = 200*60; % minutes, 2000 is more that 1 enceladus full revolution

% perturbingBodyNaifIDs = [-2, -3, -4, 606, 605, 10, 603, 604, 601, 607, 608, 602, 5];
perturbingBodyNaifIDs = [];

% pars.EncPlotSize = 1; %true enceladus size, usually too tiny
pars.EncPlotSize = 50;

%% ========================================================================
%  1. DEFINE INPUTS & INITIALIZE
%  ========================================================================

% --- Flyby and Central Body Parameters ---
pars.INPUTS.idCentral = 6;      % Central Body: Saturn
pars.INPUTS.idMoon    = 1;      % Enceladus flyby
pars.INPUTS.V_inf     = 4;    % Hyperbolic excess velocity at Enceladus [km/s]
pars.INPUTS.Flyby.min_h = 25;   % Minimum flyby altitude [km]

%%% TABLE 10 DEIMOS REPORRT %%%%

% % Starting Crank Targeting Manoeuvre
% nodein  = [4, deg2rad(8.6918), deg2rad(0)];
% nodeout = [4, deg2rad(8.6918), deg2rad(-86.9406)];
% 
% % Partial-COT 1 (O/I) - 1st Flyby
% nodein  = [4, deg2rad(8.6918), deg2rad(-86.9406)];
% nodeout = [4, deg2rad(8.6918), deg2rad(-88.1610)];
% 
% % 2nd Flyby
% nodein  = [4, deg2rad(8.6918), deg2rad(-88.1610)];
% nodeout = [4, deg2rad(8.6918), deg2rad(-89.3871)];
% 
% % 3rd Flyby
% nodein  = [4, deg2rad(8.6918), deg2rad(-89.3871)];
% nodeout = [4, deg2rad(8.6918), deg2rad(-90.6075)];
% 
% % 4th Flyby
% nodein  = [4, deg2rad(8.6918), deg2rad(-90.6075)];
% nodeout = [4, deg2rad(8.6918), deg2rad(-91.8337)];
% 
% % 5th Flyby
% nodein  = [4, deg2rad(8.6918), deg2rad(-91.8337)];
% nodeout = [4, deg2rad(8.6918), deg2rad(-93.0598)];
% 
% % Inclination Change Manoeuvre
% nodein  = [4, deg2rad(8.6918), deg2rad(-93.0598)];
% nodeout = [4, deg2rad(8.6918), deg2rad(180.0)];
% 
% % Inbound - Outbound Pseudo-Resonant Transfer
% nodein  = [4, deg2rad(8.6918), deg2rad(180.0)];
% nodeout = [4, deg2rad(8.5886), deg2rad(180.0)];
% 
% % Re-Enter 7:1 Resonance
% nodein  = [4, deg2rad(8.5886), deg2rad(0)];
% nodeout = [4, deg2rad(8.6918), deg2rad(0)];
% 
% % Starting Crank Targeting Manoeuvre
% nodein  = [4, deg2rad(8.6918), deg2rad(0)];
% nodeout = [4, deg2rad(8.6918), deg2rad(86.9406)];
% 
% % Partial-COT 2 (O/I) - 1st Flyby
% nodein  = [4, deg2rad(8.6918), deg2rad(86.9406)];
% nodeout = [4, deg2rad(8.6918), deg2rad(88.1610)];
% 
% % 2nd Flyby
% nodein  = [4, deg2rad(8.6918), deg2rad(88.1610)];
% nodeout = [4, deg2rad(8.6918), deg2rad(89.3871)];
% 
% % 3rd Flyby
% nodein  = [4, deg2rad(8.6918), deg2rad(89.3871)];
% nodeout = [4, deg2rad(8.6918), deg2rad(90.6075)];
% 
% % 4th Flyby
% nodein  = [4, deg2rad(8.6918), deg2rad(90.6075)];
% nodeout = [4, deg2rad(8.6918), deg2rad(91.8337)];
% 
% % 5th Flyby
nodein  = [4, deg2rad(8.6918), deg2rad(91.8337)];
nodeout = [4, deg2rad(8.6918), deg2rad(93.0598)];


% --- SPICE and Time Parameters ---
kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
% Load kernels in the main session for pre-loop calculations
loadSpiceKernels(kernels); 

spiceParam.frame    = 'J2000';
spiceParam.abcorr   = 'NONE';
spiceParam.observer = '699'; % Saturn Center NAIF ID

% --- EPOCH SWEEP PARAMETERS ---
epoch_start = date2mjd2000([2030 1 1 0 0 0]);
epoch_end = date2mjd2000([2033 1 1 0 0 0]);
num_epochs = 100;
epoch_array = linspace(epoch_start, epoch_end, num_epochs);

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
pars.Moon.Period = 2*pi*sqrt(pars.Moon.OrbRad^3/pars.Moon.mu);
pars.Moon.HillSph = pars.Moon.OrbRad*( pars.Moon.mu/(3*(pars.Moon.mu + pars.Planet.mu)))^(1/3);

% --- Parameters required for Linked Conic Calculation ---
pars.GroundTr.npoints = 300; 

% Calculate SOI radius for reference (still needed for plotting normalization)
m_enceladus = 1.08e20;
m_saturn = getAstroConstants('Saturn', 'Mass');
a_enceladus = 238020;
r_soi_enceladus = a_enceladus * (m_enceladus / m_saturn)^(2/5);

pars.INPUTS.Flyby.hMapping = 300;

pars.INPUTS.NAIFCentral = 699;
pars.INPUTS.NAIFMoon = 602;

%% ========================================================================
%  2. CALCULATE ENCELADUS AND TITAN ORBITS (ONCE)
%  ========================================================================

% Use a reference epoch for calculating moon orbits
reference_epoch = epoch_array(1);
pars.INPUTS.epoch0 = reference_epoch;

duration_sec = pars.GroundTr.t_prop * 60;
time_steps = pars.GroundTr.npoints;
time_vector_fwd = linspace(eps, duration_sec, time_steps)';

time_vector_fwdP = linspace(eps, duration_sec*4, time_steps)';

% --- Calculate the position history of Enceladus ---
r_enc_history = zeros(time_steps, 3);
v_enc_history = zeros(time_steps, 3);
for i = 1:time_steps
    current_mjd = reference_epoch + time_vector_fwd(i) / 86400;
    [r_enc, v_enc] = EphSS_car_spice2(602, current_mjd, true, spiceParam);
    r_enc_history(i, :) = r_enc;
    v_enc_history(i, :) = v_enc;
end

% --- Calculate the position history of Titan ---
r_titan_history = zeros(time_steps, 3);
v_titan_history = zeros(time_steps, 3);
for i = 1:time_steps
    current_mjd = reference_epoch + time_vector_fwdP(i) / 86400;
    [r_titan, v_titan] = EphSS_car_spice2(606, current_mjd, true, spiceParam);
    r_titan_history(i, :) = r_titan;
    v_titan_history(i, :) = v_titan;
end

fprintf('Ephemeris calculation completed for Enceladus and Titan.\n');

%% ========================================================================
%  3. LOOP THROUGH EPOCHS AND CALCULATE TRAJECTORIES (PARALLEL)
%  ========================================================================

% --- Start Parallel Pool and Initialize Workers ---
if isempty(gcp('nocreate'))
    parpool; % Start a parallel pool if one is not running
end
spmd % Execute on each worker
    loadSpiceKernels(kernels);
end

% [Previous code remains the same until the parfor loop section]

% Pre-allocate storage for all trajectories and Keplerian elements
all_trajectories_LC = cell(num_epochs, 1);
all_times = cell(num_epochs, 1);
all_kepga_in = zeros(num_epochs, 6);   % Keplerian elements for incoming trajectory
all_kepga_out = zeros(num_epochs, 6);  % Keplerian elements for outgoing trajectory

fprintf('Starting parallel epoch sweep from %s to %s\n', ...
    mjd20002date(epoch_start), mjd20002date(epoch_end));

parfor epoch_idx = 1:num_epochs
    % Note: fprintf progress indicators do not work as expected inside a parfor loop.
    
    % Each worker gets its own copy of 'pars' and modifies it locally.
    local_pars = pars;
    
    current_epoch = epoch_array(epoch_idx);
    local_pars.INPUTS.epoch0 = current_epoch;
    
    %% Calculate Initial Pericenter State for Propagation
    [r_enceladus_at_flyby, v_enceladus_at_flyby] = EphSS_car_spice2(602, current_epoch, true, spiceParam);
    [vvinfin, r0_sc_in, v0_sc_in, ~, kepga_in] = vinfAlphaCrank_to_VinfCARTClaudio(nodein(1), nodein(2), nodein(3), current_epoch, local_pars);
    [vvinfout, r0_sc_out, v0_sc_out, vvga, kepga_out] = vinfAlphaCrank_to_VinfCARTClaudio(nodeout(1), nodeout(2), nodeout(3), current_epoch, local_pars);
    
    rp_flyby  = local_pars.INPUTS.Flyby.min_h + local_pars.Moon.EquRad;
    e_fly     = 1 + ((rp_flyby*local_pars.INPUTS.V_inf^2)/local_pars.Moon.mu);
    delta_max = 2*asin(1/e_fly);
    % Storing this in local_pars, though it is not used later in this iteration
    local_pars.delta_max = delta_max; 
    
    Rm = buildRm(r0_sc_out,vvga);
    
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
    
    %% Propagate the Linked Conics Trajectory
    time_vector_bwd = linspace(0, -duration_sec, time_steps)';
    
    % --- Propagate backward from the incoming node ---
    [timeLC_bwd, stateLC_bwd, ~, ~, ~] = propagateKeplerODE2(r0_sc_in', v0_sc_in', time_vector_bwd, mu_central_body, []);
    
    % --- Propagate forward from the outgoing node ---
    [timeLC_fwd, stateLC_fwd, ~, ~, ~] = propagateKeplerODE2(r0_sc_out', v0_sc_out', time_vector_fwd, mu_central_body, []);
    
    % Flip the backward propagation results
    stateLC_bwd_flipped = flipud(stateLC_bwd);
    timeLC_bwd_flipped = flipud(timeLC_bwd);
    
    % Merge the two branches
    full_time_out_LC = [timeLC_bwd_flipped; timeLC_fwd(1:end)];
    full_state_out_LC = [stateLC_bwd_flipped; stateLC_fwd(1:end, :)];
    
    % Store the results for this iteration
    all_trajectories_LC{epoch_idx} = full_state_out_LC;
    all_times{epoch_idx} = full_time_out_LC;
    all_kepga_in(epoch_idx, :) = kepga_in;
    all_kepga_out(epoch_idx, :) = kepga_out;
end

fprintf('Epoch sweep completed!\n');

%% ========================================================================
%  4. CREATE INTERACTIVE VISUALIZATION WITH EPOCH SLIDER
%  ========================================================================

liveplot_epoch_sweep(all_trajectories_LC, all_times, r_enc_history, r_titan_history, epoch_array, pars, all_kepga_out);


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