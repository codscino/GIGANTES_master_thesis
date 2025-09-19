clc;
clear all;
close all;

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

% --- Define Perturbations (Saturn + Enceladus ONLY) ---
perturbingBodyNaifIDs = [602]; % enceladus
specialPerturbationIDs = []; % No J2, J3, etc. for this run

% --- Load Gravitational Parameters ---
mu_central_body = getAstroConstants('Saturn', 'Mu');
[~, mu_enceladus, ~, ~] = satMoonsConstants(1); % Enceladus
mu_TBs = mu_enceladus;

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



%% ========================================================================
%  2. PROPAGATE THE TRAJECTORY
%  ========================================================================
% --- Set Propagation Duration ---
propagation_duration_days = 100;
% propagation_duration_days = 0.001;
duration_sec = propagation_duration_days * 86400;
time_vector = linspace(0, duration_sec, 1000)'; % steps


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

%% ========================================================================
%  4. VISUALIZATION
%  ========================================================================

% --- Calculate SOI of Enceladus ---
m_enceladus = 1.08e20;
m_saturn = getAstroConstants('Saturn', 'Mass');
a_enceladus = 238020;
r_soi_enceladus = a_enceladus * (m_enceladus / m_saturn)^(2/5);

% Normalize the distance for the x-axis
distance_in_soi_units = distance_from_enceladus_history / r_soi_enceladus;

% --- Calculate the final stabilized semi-major axis ---
final_a_stabilized = mean(semi_major_axis_history(round(0.9*end):end));

% ************************************
% *** find the stabilization point ***
% ************************************
tolerance_km = 0.5; % Define the tolerance
% Find the absolute difference between each point and the final value
absolute_difference_km = abs(semi_major_axis_history - final_a_stabilized);

% Search backwards to find the last point that is OUTSIDE the tolerance
stabilization_index = -1; % Default value if all points are within tolerance
for i = length(absolute_difference_km):-1:1
    if absolute_difference_km(i) > tolerance_km
        % This is the last point to violate the condition. The stable
        % region begins at the next index.
        stabilization_index = i + 1;
        break; % Exit the loop once we've found it
    end
end

% If the whole trajectory was stable, default to the first point
if stabilization_index == -1
    stabilization_index = 1;
end
% Prevent index from going out of bounds
if stabilization_index > length(distance_in_soi_units)
    stabilization_index = length(distance_in_soi_units);
end

% Get the distance in SOI units at which stabilization occurs
stabilization_distance_soi = distance_in_soi_units(stabilization_index);
fprintf('Semi-major axis stabilizes within %.1f km at %.2f SOI radii.\n', ...
        tolerance_km, stabilization_distance_soi);
% *************************************************************************

figure('Name', 'Semi-Major Axis Evolution vs. Normalized Distance', 'Position', [100, 100, 1000, 600]);
plot(distance_in_soi_units, semi_major_axis_history / 1e6, 'b-', 'LineWidth', 2);
hold on;

% Add the vertical line for the edge of the SOI
xline(1, 'k--', 'LineWidth', 1.5, 'Label', 'Edge of SOI (x=1)');

% Add the vertical line of the tolerance met
xline(stabilization_distance_soi, 'g--', 'LineWidth', 1.5, ...
    'Label', sprintf('Stabilized (±%.1f km) at %.1f SOI', tolerance_km, stabilization_distance_soi));

% Add the horizontal line for the final value
if ~isnan(final_a_stabilized)
    yline(final_a_stabilized / 1e6, 'r--', 'LineWidth', 1.5, ...
        'Label', sprintf('Final "a" ≈ %.2f M km', final_a_stabilized/1e6));
end

title('Evolution of Saturn-Centric Semi-Major Axis vs. Distance from Enceladus');
xlabel('Distance from Enceladus (in number of SOI radii)');
ylabel('Osculating Semi-Major Axis (millions of km)');
legend('Evolving Semi-Major Axis', 'Location', 'Best');
grid on;
set(gca, 'FontSize', 12);

fprintf('Plot generated. Enceladus SOI is approx %.0f km.\n', r_soi_enceladus);

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