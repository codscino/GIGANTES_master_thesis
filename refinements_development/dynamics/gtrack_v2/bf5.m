clc;
clear all;
close all;
% script for the comparison of LC and nbody groudntracks both propagated
% from the pericentre
%% ========================================================================
%  1. DEFINE INPUTS & INITIALIZE
%  ========================================================================

% --- Flyby and Central Body Parameters ---
pars.INPUTS.idCentral = 6;      % Central Body: Saturn
pars.INPUTS.idMoon    = 1;      % Enceladus flyby
pars.INPUTS.V_inf     = 4.0;    % Hyperbolic excess velocity at Enceladus [km/s]
pars.INPUTS.Flyby.min_h = 25;   % Minimum flyby altitude [km]

% --- Define Incoming and Outgoing Asymptotes (Nodes) ---
nodein = [4, deg2rad(8.6918), deg2rad(-86.9406)];
nodeout = [4, deg2rad(8.6918), deg2rad(-88.1610)];

% --- SPICE and Time Parameters ---
kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels); 

spiceParam.frame    = 'J2000';
spiceParam.abcorr   = 'NONE';
spiceParam.observer = '699'; % Saturn Center NAIF ID
pars.INPUTS.epoch0  = date2mjd2000([2035, 1, 1, 0, 0, 0]);

% --- Define Perturbations (Saturn + Enceladus ONLY) ---
% perturbingBodyNaifIDs = [602]; % Enceladus NAIF ID
perturbingBodyNaifIDs = [-2, -3, -4, 606, 605, 10, 603, 604, 601, 607, 608, 602, 5];
% perturbingBodyNaifIDs = [-2, -3, -4, 606, 605, 10, 603, 604, 601, 607, 608];

actualBodyNaifIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs >= 0);
specialPerturbationIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs < 0);

% --- Load Gravitational Parameters ---
mu_central_body = getAstroConstants('Saturn', 'Mu');
[~, mu_enceladus, R_enceladus, ~] = satMoonsConstants(1); % Enceladus
% mu_TBs = mu_enceladus;

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
pars.GroundTr.t_prop       = 5;    % min

pars.INPUTS.Flyby.hMapping = 300;

% propagate to 64 SOI (stable sma)
m_enceladus = 1.08e20;
m_saturn = getAstroConstants('Saturn', 'Mass');
a_enceladus = 238020;
r_soi_enceladus = a_enceladus * (m_enceladus / m_saturn)^(2/5);
% pars.INPUTS.Flyby.hMapping = r_soi_enceladus*64;

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

b1 = -r0_sc_out./norm(r0_sc_out);
b3 = cross(r0_sc_out, vvga)./norm(cross(r0_sc_out, vvga));
b2 = cross(b3,b1);
Rm = [ b1' b2' b3' ]';

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
%  3. PROPAGATE THE TRAJECTORIES
%  ========================================================================

propagation_duration_mins = pars.GroundTr.t_prop;
duration_sec = propagation_duration_mins * 60;
time_steps = pars.GroundTr.npoints;
time_vector_fwd = linspace(0, duration_sec, time_steps)';
time_vector_bwd = linspace(0, -duration_sec, time_steps)';

ephem_handle = @(t_sec) get_body_positions_wrapper(t_sec, pars.INPUTS.epoch0, actualBodyNaifIDs, spiceParam);
[~, state_out_fwd] = propagateNBodyODE2(rrp_saturn_centric, vvp_saturn_centric, time_vector_fwd, mu_central_body, mu_TBs, ephem_handle, specialPerturbationIDs);
[time_out_bwd, state_out_bwd] = propagateNBodyODE2(rrp_saturn_centric, vvp_saturn_centric, time_vector_bwd, mu_central_body, mu_TBs, ephem_handle, specialPerturbationIDs);

time_out_bwd = flipud(time_out_bwd);
state_out_bwd = flipud(state_out_bwd);
full_time_out = [time_out_bwd; time_vector_fwd(2:end)];
full_state_out = [state_out_bwd; state_out_fwd(2:end, :)];
num_steps = length(full_time_out);


%% ========================================================================
%  4. PREPARE DATA FOR PLOTTING
%  ========================================================================

% --- Calculate the full position history of Enceladus ---
r_enc_history = zeros(num_steps, 3);
v_enc_history = zeros(num_steps, 3);

for i = 1:num_steps
    current_mjd = pars.INPUTS.epoch0 + full_time_out(i) / 86400;
    [r_enc, v_enc] = EphSS_car_spice2(602, current_mjd, true, spiceParam);
    r_enc_history(i, :) = r_enc;
    v_enc_history(i, :) = v_enc;
end

% --- Compute N-Body ground track ---

lats = zeros(num_steps, 1);
lon_rad_temp = zeros(num_steps, 1);
alts = zeros(num_steps, 1);
r_in_flyby_frame = zeros(num_steps, 3);

for i = 1:num_steps
    current_mjd = pars.INPUTS.epoch0  + full_time_out(i) / 86400;
    sc_pos_icrf = full_state_out(i, 1:3); 
    r_enc = r_enc_history(i, :);
    v_enc = v_enc_history(i, :);

    % --- TRANSFORMATION USING a time-varying Rm ---
    sc_pos_relative_icrf = sc_pos_icrf - r_enc;
    Rm = buildRm(r_enc, v_enc);
    r_in_flyby_frame(i, :) = (Rm * sc_pos_relative_icrf')';
    r_in_flyby_frame_unit = r_in_flyby_frame(i, :) / norm(r_in_flyby_frame(i, :));
  

    % % --- TRANSFORMATION USING Rm costant approximation ---
    % 
    % % 1. Translate: Get S/C position relative to Enceladus(still ICRF)
    % sc_pos_relative_icrf = sc_pos_icrf - r_enc;
    % 
    % % This transforms the vector into the static flyby geometry frame.
    % r_in_flyby_frame = (Rm * sc_pos_relative_icrf')';
    % 
    % 3. Calculate "latitude" and "longitude" based on this static frame.
    % and nomralize the vector to get the angles correctly.
    % r_in_flyby_frame_unit = r_in_flyby_frame / norm(r_in_flyby_frame);

    lon_rad_temp(i) = atan2(r_in_flyby_frame_unit(2), r_in_flyby_frame_unit(1));
    lats(i) = asin(r_in_flyby_frame_unit(3));
    
    % Altitude is a scalar distance, so it is calculated from the relative
    % vector before rotation and is unaffected by the approximation.
    alts(i) = norm(sc_pos_relative_icrf) - pars.Moon.EquRad;
end

 % --> make it comparable with Campagnola --> shift to east-longitude
 long   = wrapTo2Pi(lon_rad_temp);
 longs = wrapTo360(long);

 % convert to degrees
lats_deg = rad2deg(lats);
longs_deg = rad2deg(longs);

% finds the the index of the periapsis by finding the minimum altitude
[~, idx_periapsis] = min(alts); % if time_steps = 30e3; then idx_periapsis=30000

% --- Compute Linked Conic reference ground track ---
[Flyby_LinkedConic] = Flyby_BuildUp_claudio(nodein, nodeout, pars); %spice

% [Flyby_LinkedConic] = Flyby_BuildUp(nodein, nodeout, pars); %no spice


%% ========================================================================
%  5. PLOT TRAJECTORIES AND GROUND TRACKS
%  ========================================================================
 
% --- Create a new, wider figure with a white background ---
figure('Name', 'Trajectory and Ground Track Analysis', 'Color', 'w', 'Position', [100 100 1400 600]);

% --- LEFT SUBPLOT: 3D Static Trajectory Plot ---
ax1 = subplot(1, 2, 1);
hold(ax1, 'on');

% Plot the full trajectories
plot3(ax1, full_state_out(:,1), full_state_out(:,2), full_state_out(:,3), 'b-', 'LineWidth', 2, 'DisplayName', 'Spacecraft');
plot3(ax1, r_enc_history(:,1), r_enc_history(:,2), r_enc_history(:,3), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Enceladus');

% --- Define colors and styles from Plot_Flyby_GT for consistency ---
color_exit = [0.9290 0.6940 0.1250]; % Gold/Orange for exit
marker_size_3d = 6;

% Plot start, end, and pericenter points with matching styles
plot3(ax1, full_state_out(1,1), full_state_out(1,2), full_state_out(1,3), ...
    'o', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black', 'MarkerSize', marker_size_3d, ...
    'DisplayName', 'Flyby Entry Point');

plot3(ax1, full_state_out(end,1), full_state_out(end,2), full_state_out(end,3), ...
    'o', 'MarkerFaceColor', color_exit, 'MarkerEdgeColor', 'black', 'MarkerSize', marker_size_3d, ...
    'DisplayName', 'Flyby Exit Point');

plot3(ax1, rrp_saturn_centric(1), rrp_saturn_centric(2), rrp_saturn_centric(3), ...
    'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black', 'MarkerSize', marker_size_3d, ...
    'DisplayName', 'Flyby Periapsis');

% Formatting
title(ax1, 'Spacecraft and Enceladus Trajectory (ICRF)');
xlabel(ax1, 'X [km]');
ylabel(ax1, 'Y [km]');
zlabel(ax1, 'Z [km]');
legend(ax1, 'show', 'Location', 'northeast');
grid(ax1, 'on');
axis(ax1, 'equal');
view(ax1, 3);

% --- RIGHT SUBPLOT: Ground Track Comparison ---
ax2 = subplot(1, 2, 2);
hold(ax2, 'on');

% --- Filter data for mapping altitude ---
mapping_altitude_km = pars.INPUTS.Flyby.hMapping;
valid_indices = alts <= mapping_altitude_km;
NBody_Mapping_Data.lats = lats(valid_indices);
NBody_Mapping_Data.longs = longs(valid_indices);
NBody_Mapping_Data.rp_lat = lats_deg(idx_periapsis);
NBody_Mapping_Data.rp_long = longs_deg(idx_periapsis);

% Plot the moon texture map
plotTextureLatLong(pars.INPUTS.idMoon, pars.INPUTS.idCentral, 1); 
% xticks(0:15:360);
% xticklabels({'0°E','15°E','30°E','45°E','60°E','75°E','90°E','105°E','120°E','135°E','150°E','165°E',...
%             '180°E','195°E','210°E','225°E','240°E','255°E','270°E','285°E','300°E','315°E','330°E','345°E','360°E'});

axis(ax2, 'normal');
daspect(ax2, [1 1 1]);  % Equal data aspect ratio
grid(ax2, 'on');
title(ax2, ['Ground Track Comparison (Altitude \leq ' num2str(mapping_altitude_km) ' km)']);

% Plot Linked Conic Ground Track (Blue)
if ~isempty(Flyby_LinkedConic)
     Plot_Flyby_GTClaudio(Flyby_LinkedConic, [0, 0.4470, 0.7410], '--', 4);
end

% Plot N-Body Ground Track (Orange)
if ~isempty(NBody_Mapping_Data.lats)
    Plot_Flyby_GTClaudio(NBody_Mapping_Data, [0.8500, 0.3250, 0.0980 0.7]); 
end



% Reset the legend to only show the main track types, as Plot_Flyby_GT creates its own detailed legends.
h = zeros(2, 1);
h(1) = plot(ax2, NaN,NaN,'-', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2);
h(2) = plot(ax2, NaN,NaN,'-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2);
legend(h, 'N-Body Propagation', 'Linked Conic Approx.', 'Location', 'best');

hold(ax1, 'off');
hold(ax2, 'off');


%% ========================================================================
%  6. PLOT TRAJECTORIES IN ENCELADUS-CENTERED STATIC FRAME
%  ========================================================================

% Get the total number of time steps from the output state matrix
num_steps = size(full_state_out, 1);

% Pre-allocate memory for the transformed position vectors for efficiency
r_sc_in_flyby_frame = zeros(num_steps, 3);

% Loop through each time step of the propagated trajectory
for i = 1:num_steps
    % --- Step 1: Translation ---
    % Get the spacecraft's position relative to Enceladus by subtracting
    % Enceladus's position from the spacecraft's position.
    % The result is still in the ICRF (J2000) frame.
    sc_pos_relative_icrf = full_state_out(i, 1:3) - r_enc_history(i, :);

    % --- Step 2: Rotation ---
    % Rotate the relative position vector from the ICRF frame into the
    % static flyby frame defined by the constant matrix Rm.
    % Note: The position vector must be transposed to a column vector for matrix multiplication.
    r_sc_in_flyby_frame(i, :) = (Rm * sc_pos_relative_icrf')';
end

r_sc_in_flyby_frameLC = Flyby_LinkedConic.fly_States(:,1:3);


% --- Create a new figure for the 3D plot ---
figure('Name', 'Enceladus Flyby in Static Frame', 'Color', 'w');
hold on; % Hold the plot to draw multiple items

% --- Plot Enceladus as a sphere ---
% Define the radius and a light blue color for the moon
radius = pars.Moon.EquRad;
light_blue_color = [0.678, 0.847, 0.902]; % light blue

% Generate the coordinates for a sphere
[x_sphere, y_sphere, z_sphere] = sphere(50); 

% Draw the sphere, scaling it by the moon's radius
surf(x_sphere*radius, y_sphere*radius, z_sphere*radius, ...
     'FaceColor', light_blue_color, ...
     'EdgeColor', 'none', ...
     'DisplayName', 'Enceladus');
alpha(0.7); % Make the sphere slightly transparent to see trajectories behind it

% --- Plot the Linked Conics Trajectory ---
plot3(r_sc_in_flyby_frameLC(:,1), r_sc_in_flyby_frameLC(:,2), r_sc_in_flyby_frameLC(:,3), ...
      'Color', [0, 0.4470, 0.7410], ... % Blue color
      'LineStyle', '--', ...             % Dashed line
      'LineWidth', 2, ...
      'DisplayName', 'Linked Conic Trajectory');

% --- Plot the N-Body Trajectory ---
plot3(r_sc_in_flyby_frame(:,1), r_sc_in_flyby_frame(:,2), r_sc_in_flyby_frame(:,3), ...
      'Color', [0.8500, 0.3250, 0.0980 0.6], ... % Orange color
      'LineWidth', 4, ...
      'DisplayName', 'N-Body Trajectory');

% --- Add plot formatting and labels ---
title('Flyby Trajectory in Enceladus-Centered Static Frame');
xlabel('Flyby Frame X [km]');
ylabel('Flyby Frame Y [km]');
zlabel('Flyby Frame Z [km]');
legend('show', 'Location', 'best'); % Display the legend
grid on;       % Add a grid
axis equal;    % Ensure the sphere looks like a sphere
view(3);       % Set the default view to 3D

hold off; 


%% ========================================================================
%  7. DEVIATION ANALYSIS: LINKED CONICS VS. N-BODY
%  ========================================================================


% --- Calculate the Euclidean distance (norm of the difference) between the
% --- N-body and Linked Conic position vectors at each time step.
% --- Both trajectories are in the same Enceladus-centered static frame.
deviation_km = vecnorm(r_sc_in_flyby_frame - r_sc_in_flyby_frameLC, 2, 2);

% --- Create a new figure for the deviation plots ---
figure('Name', 'N-Body vs. Linked Conic Deviation Analysis', 'Color', 'w', 'Position', [200 200 1200 500]);

% --- SUBPLOT 1: Deviation as a function of Time ---
ax1 = subplot(1, 2, 1);
hold(ax1, 'on');

% Convert time vector from seconds to minutes for better readability
time_in_minutes = full_time_out / 60;

% Plot the full deviation over time
plot(ax1, time_in_minutes, deviation_km, 'b-', 'LineWidth', 2, 'DisplayName', 'Position Deviation');

% Highlight the portion of the trajectory below the mapping altitude
% This corresponds to the 'valid_indices' calculated earlier
plot(ax1, time_in_minutes(valid_indices), deviation_km(valid_indices), ...
    'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2.5, 'DisplayName', ['Deviation < ' num2str(mapping_altitude_km) ' km Alt']);

% Add a vertical line to mark the time of periapsis (t=0)
xline(ax1, 0, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Periapsis (t=0)', 'Label', 'Periapsis');

% Formatting
title(ax1, 'Trajectory Deviation vs. Time');
xlabel(ax1, 'Time from Periapsis [minutes]');
ylabel(ax1, 'Position Deviation [km]');
legend(ax1, 'show', 'Location', 'northwest');
grid(ax1, 'on');
hold(ax1, 'off');


% --- SUBPLOT 2: Deviation as a function of Altitude ---
ax2 = subplot(1, 2, 2);

% Plot deviation against altitude 
plot(ax2, alts, deviation_km, 'k.', 'MarkerSize', 0.1, 'DisplayName', 'All Data Points'); 
hold(ax2, 'on');

% red markers under 300km
scatter(ax2, alts(valid_indices), deviation_km(valid_indices), 5, 'r', 'filled', 'DisplayName', ['Altitude < ' num2str(mapping_altitude_km) ' km']); 

% Formatting
title(ax2, 'Trajectory Deviation vs. Altitude');
xlabel(ax2, 'Altitude Above Enceladus [km]');
ylabel(ax2, 'Position Deviation [km]');
legend(ax2, 'show', 'Location', 'best');
grid(ax2, 'on');

% Set x-axis to log scale to better visualize the behavior at both
% close and far distances from the moon.
set(ax2, 'XScale', 'log');
hold(ax2, 'off');



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

