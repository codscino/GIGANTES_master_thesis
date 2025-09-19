%%% draw groundtracks with nbody propagator %%%
clc;
clear all;
close all;


%% soi earth
m_earth = getAstroConstants('Earth', 'mass');
m_sun = getAstroConstants('Sun', 'mass');
au = getAstroConstants('AU');

r_soi = au*(m_earth/m_sun)^(2/5);
r_soi_au = r_soi/au;

n_sois = 0.4/r_soi_au;

%% ========================================================================
%  1. DEFINE INPUTS & INITIALIZE
%  ========================================================================

% --- Flyby and Central Body Parameters ---
pars.INPUTS.idCentral = 6;      % Central Body: Saturn
pars.INPUTS.idMoon    = 1;      % Enceladus flyby
pars.INPUTS.V_inf     = 4.0;    % Hyperbolic excess velocity at Enceladus [km/s]
pars.INPUTS.k         = 0;      % Crank angle [rad] (fixed)
pars.INPUTS.alpha     = 0.15;   % Pump angle [rad] (fixed)

pars.INPUTS.Flyby.min_h    = 25;

% --- SPICE and Time Parameters ---
kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels);

spiceParam.frame    = 'J2000';
spiceParam.abcorr   = 'NONE';
spiceParam.observer = '699'; % Saturn Center NAIF ID

t0_mjd = date2mjd2000([2035, 1, 1, 0, 0, 0]);

all_perturber_ids = [602, -2, -4, 606, 605, 10, 603, 604];

% --- Filter into actual bodies and special perturbations ---
perturbingBodyNaifIDs = all_perturber_ids(all_perturber_ids >= 0);
specialPerturbationIDs = all_perturber_ids(all_perturber_ids < 0);

% --- Load Gravitational Parameter for the Central Body ---
mu_central_body = getAstroConstants('Saturn', 'Mu');

% --- Dynamically build the list of gravitational parameters (mu) for all perturbing bodies ---
mu_TBs = zeros(1, length(perturbingBodyNaifIDs));
for i = 1:length(perturbingBodyNaifIDs)
    id = perturbingBodyNaifIDs(i);
    switch id
        case 10,  mu_TBs(i) = getAstroConstants('Sun', 'Mu');
        case 5,   mu_TBs(i) = getAstroConstants('Jupiter', 'Mu');
        case 601, [~, mu, ~, ~] = satMoonsConstants(0); mu_TBs(i) = mu; % Mimas
        case 602, [~, mu, ~, ~] = satMoonsConstants(1); mu_TBs(i) = mu; % Enceladus
        case 603, [~, mu, ~, ~] = satMoonsConstants(2); mu_TBs(i) = mu; % Tethys
        case 604, [~, mu, ~, ~] = satMoonsConstants(3); mu_TBs(i) = mu; % Dione
        case 605, [~, mu, ~, ~] = satMoonsConstants(4); mu_TBs(i) = mu; % Rhea
        case 606, [~, mu, ~, ~] = satMoonsConstants(5); mu_TBs(i) = mu; % Titan
        case 607, mu_TBs(i) = 0.374;   % Hyperion
        case 608, mu_TBs(i) = 120.4;   % Iapetus
        otherwise
            error('Unknown NAIF ID %d for mu calculation.', id);
    end
end


% --- Pericentre starting point ---
[~, r0_sc, v0_sc, ~] = vinfAlphaCrank_to_VinfCARTClaudio(pars.INPUTS.V_inf, ...
    pars.INPUTS.alpha, pars.INPUTS.k, t0_mjd);


% otherwise propagator divides by zero because sc and encelaus are
% in the same position
r_enc_unit = r0_sc / norm(r0_sc);
r_offset_direction = [-r_enc_unit(2); r_enc_unit(1); 0]; 
rp_flyby  = pars.INPUTS.Flyby.min_h + 252;   
r_offset = r_offset_direction * rp_flyby;
r0_sc = r0_sc + r_offset';

% Ensure vectors are columns for the propagator
r0_sc = r0_sc(:);
v0_sc = v0_sc(:);

% --- Set Propagation Duration ---
propagation_duration_days = 30; %0.001;
duration_sec = propagation_duration_days * 86400;
time_vector = linspace(0, duration_sec, 500)'; % steps

%% ========================================================================
%  2. PROPAGATE THE TRAJECTORY
%  ========================================================================

% --- Setup Ephemeris Handle using the new wrapper function ---
% This handle correctly interfaces with propagateNBodyODE2.
ephem_handle = @(t_sec) get_body_positions_wrapper(t_sec, t0_mjd, ...
    perturbingBodyNaifIDs, spiceParam);

% --- Run the N-Body Propagator ---
[time_out, state_out] = propagateNBodyODE2(r0_sc, v0_sc, time_vector, ...
    mu_central_body, mu_TBs, ephem_handle, specialPerturbationIDs);


%% ========================================================================
%  3. POST-PROCESSING: CALCULATE OSCULATING SEMI-MAJOR AXIS
%  ========================================================================
num_steps = length(time_out);
semi_major_axis_history = zeros(num_steps, 1);
distance_from_enceladus_history = zeros(num_steps, 1);

for i = 1:num_steps
    current_sc_state = state_out(i, :);

    current_mjd = t0_mjd + time_out(i) / 86400;
    
    % Get Enceladus's position relative to Saturn
    [r_enceladus_now, ~] = EphSS_car_spice2(602,...
        current_mjd, true, spiceParam);

    % --- Calculate Semi-Major Axis (relative to Saturn) ---
    kep = car2kep(current_sc_state, mu_central_body);
    semi_major_axis_history(i) = kep(1);

     % --- Calculate Distance from Enceladus ---
    r_sc_now = current_sc_state(1:3)'; % Get S/C position as column
    r_relative = r_sc_now - r_enceladus_now(:); % Vector from Enceladus to S/C
    distance_from_enceladus_history(i) = norm(r_relative);
end




% start from the pericentre and backpropagate with nbody
% convert cartesians to lat,lon on enceladus
% plot lat lon mercator projection



% plot the same thing but with the standard script and confront the 2
% gtracks


% repeat the same but with forward propagation


% repeat the same but with patched conics(enceladus 2 body eq till border
% of soi, then saturn becomes the main pertuber)

% repeat the same but with 3bodyrestricted(nbody with just enceladus added
% to the main saturn)


%% ========================================================================
%  HELPER FUNCTION
%  ========================================================================
function positions = get_body_positions_wrapper(t_sec, t0_mjd_start, naif_ids, sp_param)
% This function acts as a bridge between the propagator's requirements
% and the EphSS_car_spice2 function.
%
% Propagator expects: A 3xM matrix of positions.
% EphSS_car_spice2 provides: State of ONE body at a time.
%
    current_mjd = t0_mjd_start + t_sec / 86400;
    num_bodies = length(naif_ids);
    positions = zeros(3, num_bodies);
    
    for k = 1:num_bodies
        % Call the user-provided function to get the state of the k-th body
        [r, ~] = EphSS_car_spice2(naif_ids(k), current_mjd, true, sp_param);
        positions(:, k) = r; % Store the position vector as a column
    end
end