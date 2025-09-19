clc;
clear all;
close all;

% =========================================================================
% NOTE: This script assumes the following functions exist in the MATLAB path:
% - icrf2enceladus.m, enceladus2icrf.m, enceladus_orientation.m
% - Rotx.m, Rotz.m
% - buildRm.m (for initial flyby setup only)
% - createAnimatedFlybyPlot.m
% - createAnimatedEnceladusCentricPlot.m
% - And all other required simulation/utility functions...
% =========================================================================


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

perturbingBodyNaifIDs = [-2, -3, -4, 606, 605, 10, 603, 604, 601, 607, 608, 602, 5];
% perturbingBodyNaifIDs = [602];

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
pars.GroundTr.t_prop       = 10*2000;    % min, more that one enceladus full revolution


% propagate to 64 SOI (stable sma)
m_enceladus = 1.08e20;
m_saturn = getAstroConstants('Saturn', 'Mass');
a_enceladus = 238020;
r_soi_enceladus = a_enceladus * (m_enceladus / m_saturn)^(2/5);
pars.INPUTS.Flyby.hMapping = r_soi_enceladus*5064;


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

% Using buildRm here is part of the initial problem setup to define the
% flyby geometry. The frame it creates is temporary and is not used for
% the time-series analysis and plotting later in the script.
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
soi_event_func = @(t, x) soiCrossingEvent(t, x, ephem_enceladus, pars.INPUTS.Flyby.hMapping);

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
%  3.2 PROPAGATE THE LINKED CONICS TRAJECTORY (WITH EVENT DETECTION)
%  ========================================================================

% --- Propagate backward from the incoming node until the event ---
[timeLC_bwd, stateLC_bwd, ~, ~, ~] = propagateKeplerODE2(r0_sc_in', v0_sc_in', time_vector_bwd, mu_central_body, soi_event_func);

% --- Propagate forward from the outgoing node until the event ---
[timeLC_fwd, stateLC_fwd, ~, ~, ~] = propagateKeplerODE2(r0_sc_out', v0_sc_out', time_vector_fwd, mu_central_body, soi_event_func);

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
state_LC_resampled = interp1(full_time_out_LC, full_state_out_LC, full_time_out_nb, 'spline');

%% ========================================================================
%  4. PREPARE DATA FOR PLOTTING
%  ========================================================================

% --- Calculate the full position history of Enceladus ---
r_enc_history = zeros(num_steps, 3);
v_enc_history = zeros(num_steps, 3);
mjd_history = zeros(num_steps, 1);
for i = 1:num_steps
    current_mjd = pars.INPUTS.epoch0 + full_time_out_nb(i) / 86400;
    mjd_history(i) = current_mjd;
    [r_enc, v_enc] = EphSS_car_spice2(602, current_mjd, true, spiceParam);
    r_enc_history(i, :) = r_enc;
    v_enc_history(i, :) = v_enc;
end


%% ========================================================================
%  6. PLOT TRAJECTORIES IN IAU_ENCELADUS BODY-FIXED FRAME
%  ========================================================================

% Pre-allocate memory for the transformed position vectors
r_sc_in_iau_nb = zeros(num_steps, 3);
r_sc_in_iau_lc = zeros(num_steps, 3);

fprintf('Transforming trajectories to IAU_ENCELADUS frame using external function...\n');
for i = 1:num_steps
    % Get spacecraft ICRF position for both trajectories
    r_sc_icrf_nb = full_state_out_nb(i, 1:3);
    r_sc_icrf_lc = state_LC_resampled(i, 1:3);
    
    % Get Enceladus ICRF position at the current time
    r_enc_current = r_enc_history(i, :);
    
    % Get current time in MJD2000
    T_current = mjd_history(i);
    
    % --- Transform N-Body Trajectory ---
    dist_nb = norm(r_sc_icrf_nb - r_enc_current);
    u_sc_iau_nb = icrf2enceladus(r_sc_icrf_nb, r_enc_current, T_current);
    r_sc_in_iau_nb(i, :) = (u_sc_iau_nb * dist_nb)';
    
    % --- Transform Linked Conic Trajectory ---
    dist_lc = norm(r_sc_icrf_lc - r_enc_current);
    u_sc_iau_lc = icrf2enceladus(r_sc_icrf_lc, r_enc_current, T_current);
    r_sc_in_iau_lc(i, :) = (u_sc_iau_lc * dist_lc)';
end

% --- Create a new figure for the 3D plot ---
figure('Name', 'Enceladus Flyby in IAU Body-Fixed Frame', 'Color', 'w');
hold on; 

% --- Plot Enceladus as a sphere ---
radius = pars.Moon.EquRad;
light_blue_color = [0.678, 0.847, 0.902]; 
[x_sphere, y_sphere, z_sphere] = sphere(50); 
surf(x_sphere*radius, y_sphere*radius, z_sphere*radius, ...
     'FaceColor', light_blue_color, 'EdgeColor', 'none', 'DisplayName', 'Enceladus');
alpha(0.7); 

% --- Plot the Linked Conics Trajectory ---
plot3(r_sc_in_iau_lc(:,1), r_sc_in_iau_lc(:,2), r_sc_in_iau_lc(:,3), ...
      'Color', [0, 0.4470, 0.7410], 'LineStyle', '--', 'LineWidth', 2, ...
      'DisplayName', 'Linked Conic Trajectory');

% --- Plot the N-Body Trajectory ---
plot3(r_sc_in_iau_nb(:,1), r_sc_in_iau_nb(:,2), r_sc_in_iau_nb(:,3), ...
      'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 4, ...
      'DisplayName', 'N-Body Trajectory');

% --- Add plot formatting and labels ---
title('Flyby Trajectory in IAU_ENCELADUS Body-Fixed Frame');
xlabel('X_{IAU} [km]');
ylabel('Y_{IAU} [km]');
zlabel('Z_{IAU} [km]');
legend('show', 'Location', 'best'); 
grid on;       
axis equal;    
view(135, 25); 
camlight left;  
lighting gouraud;
hold off;


%% ========================================================================
%  7. DEVIATION ANALYSIS: LINKED CONICS VS. N-BODY
%  ========================================================================

% --- Calculate deviation in the IAU_ENCELADUS frame ---
deviation_km = vecnorm(r_sc_in_iau_nb - r_sc_in_iau_lc, 2, 2);

% --- Calculate Keplerian elements for both trajectories at each time step ---
sma_nb = zeros(num_steps, 1);
sma_lc = zeros(num_steps, 1);
for i = 1:num_steps
    kep_nb = car2kep(full_state_out_nb(i, :), mu_central_body);
    sma_nb(i) = kep_nb(1);
    kep_lc = car2kep(state_LC_resampled(i, :), mu_central_body);
    sma_lc(i) = kep_lc(1);
end
sma_difference_km = sma_nb - sma_lc;

% --- Create a new figure for the deviation plots ---
figure('Name', 'N-Body vs. Linked Conic Deviation Analysis', 'Color', 'w', 'Position', [100 100 1200 900]);
time_in_minutes = full_time_out_nb / 60;

% --- SUBPLOT 1: Position Deviation as a function of Time ---
ax1 = subplot(2, 2, 1);
hold(ax1, 'on');
plot(ax1, time_in_minutes, deviation_km, 'b-', 'LineWidth', 2, 'DisplayName', 'Position Deviation');
xline(ax1, 0, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Periapsis (t=0)', 'Label', 'Periapsis');
title(ax1, 'Position Deviation vs. Time');
xlabel(ax1, 'Time from Periapsis [minutes]');
ylabel(ax1, 'Position Deviation [km]');
legend(ax1, 'show', 'Location', 'northwest');
grid(ax1, 'on');
hold(ax1, 'off');

% --- SUBPLOT 2: Position Deviation as a function of Distance ---
ax2 = subplot(2, 2, 2);
hold(ax2, 'on');
% Calculate distance from Enceladus, normalized by its SOI radius
distance_from_enceladus_normalized = vecnorm(r_sc_in_iau_nb, 2, 2) / r_soi_enceladus;
plot(ax2, distance_from_enceladus_normalized, deviation_km, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2, 'DisplayName', 'Position Deviation');

% Mark periapsis
[min_dist, min_idx] = min(distance_from_enceladus_normalized);
dev_at_min_dist = deviation_km(min_idx);
plot(ax2, min_dist, dev_at_min_dist, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Periapsis');

% Add labels for Inbound and Outbound branches
text(ax2, distance_from_enceladus_normalized(1), deviation_km(1), ' Inbound', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
text(ax2, distance_from_enceladus_normalized(end), deviation_km(end), ' Outbound', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');

% Formatting
title(ax2, 'Position Deviation vs. Distance from Enceladus');
xlabel(ax2, 'Distance from Enceladus [Enceladus SOI]');
ylabel(ax2, 'Position Deviation [km]');
legend(ax2, 'show', 'Location', 'best');
grid(ax2, 'on');
hold(ax2, 'off');

% --- SUBPLOT 3: Semi-Major Axis Difference vs. Time ---
ax3 = subplot(2, 2, 3);
hold(ax3, 'on');
plot(ax3, time_in_minutes, sma_difference_km, 'g-', 'LineWidth', 2, 'DisplayName', 'SMA Difference');
xline(ax3, 0, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Periapsis (t=0)', 'Label', 'Periapsis');
title(ax3, 'Semi-Major Axis Difference vs. Time');
xlabel(ax3, 'Time from Periapsis [minutes]');
ylabel(ax3, 'SMA deviation [km]');
legend(ax3, 'show', 'Location', 'northwest');
grid(ax3, 'on');
hold(ax3, 'off');

% --- SUBPLOT 4: Semi-Major Axis Difference vs. Distance ---
ax4 = subplot(2, 2, 4);
hold(ax4, 'on');
plot(ax4, distance_from_enceladus_normalized, sma_difference_km, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2, 'DisplayName', 'SMA Difference');

% Mark periapsis
sma_diff_at_min_dist = sma_difference_km(min_idx);
plot(ax4, min_dist, sma_diff_at_min_dist, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Periapsis');

% Add labels for Inbound and Outbound branches
text(ax4, distance_from_enceladus_normalized(1), sma_difference_km(1), ' Inbound', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
text(ax4, distance_from_enceladus_normalized(end), sma_difference_km(end), ' Outbound', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');

% Formatting
title(ax4, 'Semi-Major Axis Difference vs. Distance from Enceladus');
xlabel(ax4, 'Distance from Enceladus [Enceladus SOI]');
ylabel(ax4, 'SMA deviation [km]');
legend(ax4, 'show', 'Location', 'best');
grid(ax4, 'on');
hold(ax4, 'off');

%% ========================================================================
%  8. PLOT ANGULAR MOMENTUM IN THE IAU_ENCELADUS BODY-FIXED FRAME
%  ========================================================================

% Calculate velocity in the rotating frame using finite differences
dt_vec = diff(full_time_out_nb);
v_sc_in_iau_nb = diff(r_sc_in_iau_nb) ./ dt_vec;

% Calculate specific angular momentum (h = r x v). Note we must use r from
% the midpoint of the interval to match the differentiated velocity.
h_iau_body = cross(r_sc_in_iau_nb(1:end-1,:), v_sc_in_iau_nb, 2);

% Create a new figure for the angular momentum plots
figure('Name', 'Specific Angular Momentum in IAU Frame', 'Color', 'w', 'Position', [300 300 1000 800]);
sgtitle('Specific Angular Momentum (h = r x v) in IAU_ENCELADUS Frame', 'FontWeight', 'bold');

% --- SUBPLOT 1: Hx component vs. Time ---
ax_hx = subplot(3, 1, 1);
plot(ax_hx, time_in_minutes(1:end-1), h_iau_body(:,1), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
title(ax_hx, 'X-Component'); ylabel(ax_hx, 'h_x [km^2/s]'); grid on;
xline(ax_hx, 0, 'r--');

% --- SUBPLOT 2: Hy component vs. Time ---
ax_hy = subplot(3, 1, 2);
plot(ax_hy, time_in_minutes(1:end-1), h_iau_body(:,2), 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2);
title(ax_hy, 'Y-Component'); ylabel(ax_hy, 'h_y [km^2/s]'); grid on;
xline(ax_hy, 0, 'r--');

% --- SUBPLOT 3: Hz component vs. Time ---
ax_hz = subplot(3, 1, 3);
plot(ax_hz, time_in_minutes(1:end-1), h_iau_body(:,3), 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2);
title(ax_hz, 'Z-Component'); ylabel(ax_hz, 'h_z [km^2/s]'); grid on;
xlabel(ax_hz, 'Time from Periapsis [minutes]');
xline(ax_hz, 0, 'r--');

% Link the x-axes
linkaxes([ax_hx, ax_hy, ax_hz], 'x');
xlim(ax_hx, [time_in_minutes(1), time_in_minutes(end-1)]);


%% ========================================================================
%  9. CREATE ANIMATED PLOTS
%  ========================================================================

fprintf('Preparing data for animations...\n');

% --- Prepare data for Enceladus-centric animation ---
% Your 'createAnimatedEnceladusCentricPlot' needs Saturn's position in the
% Enceladus-centric frame.
r_saturn_in_iau_nb = zeros(num_steps, 3);
for i = 1:num_steps
    % Saturn's ICRF position is the origin of the whole system
    r_saturn_icrf = [0, 0, 0]; 
    r_enc_current = r_enc_history(i, :);
    T_current = mjd_history(i);
    
    % Find the vector from Enceladus to Saturn
    dist_saturn = norm(r_saturn_icrf - r_enc_current);
    % Transform this direction into the IAU frame
    u_saturn_iau = icrf2enceladus(r_saturn_icrf, r_enc_current, T_current);
    % Scale the unit vector to get the full position vector
    r_saturn_in_iau_nb(i, :) = (u_saturn_iau * dist_saturn)';
end

% --- Call the external animation functions ---
fprintf('Starting animations...\n');
createAnimatedFlybyPlot(full_state_out_nb, state_LC_resampled, r_enc_history, full_time_out_nb, pars);
createAnimatedEnceladusCentricPlot(r_sc_in_iau_nb, r_sc_in_iau_lc, r_saturn_in_iau_nb, full_time_out_nb, pars);


%% ========================================================================
%  HELPER FUNCTIONS (Internal to this script)
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