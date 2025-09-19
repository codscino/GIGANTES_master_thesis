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
pars.INPUTS.Flyby.min_h = 25;   % Minimum flyby altitude [km]

% --- Define Incoming and Outgoing Asymptotes (Nodes) ---
% These define the direction of the hyperbolic excess velocity vector.
% [V_infinity, alpha, crank_angle]
nodein =  [pars.INPUTS.V_inf, 0.15, 0];            % [km/s, rad, rad]
nodeout = [pars.INPUTS.V_inf, 0.15, deg2rad(1)];   % [km/s, rad, rad]

% --- SPICE and Time Parameters ---
kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels); 

spiceParam.frame    = 'J2000';
spiceParam.abcorr   = 'NONE';
spiceParam.observer = '699'; % Saturn Center NAIF ID

pars.INPUTS.epoch0  = date2mjd2000([2035, 1, 1, 0, 0, 0]); % MJD2000 time of flyby

% --- Define Perturbations (Saturn + Enceladus ONLY) ---
perturbingBodyNaifIDs = [602]; % Enceladus NAIF ID

% --- Load Gravitational Parameters ---
mu_central_body = getAstroConstants('Saturn', 'Mu');
[~, mu_enceladus, R_enceladus, ~] = satMoonsConstants(1); % Enceladus
mu_TBs = mu_enceladus;


% Retrieve Saturn Parameters 
[pars.Planet.mu, pars.Planet.EquRad, pars.Planet.OrbRad, pars.Planet.hmin] = planetConstants(pars.INPUTS.idCentral); %[km3/s2],[km],[km] & [km]

% Retrieve Desired Moon Parameters 
if pars.INPUTS.idCentral == 3
    pars.Moon.OrbRad = 384748; pars.Moon.mu  = getAstroConstants('Moon','mu'); %[km],[km3/s2]
    pars.Moon.EquRad = getAstroConstants('Moon','Radius'); pars.Moon.hmin = 50;  %[km], [km]
elseif pars.INPUTS.idCentral == 5
    [pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = jupMoonsConstants(pars.INPUTS.idMoon); %[km],[km3/s2],[km] & [km]
elseif pars.INPUTS.idCentral == 6
    [pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = satMoonsConstants(pars.INPUTS.idMoon); %[km],[km3/s2],[km] & [km]
else
    [pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = uranusMoonsConstants(pars.INPUTS.idMoon); %[km],[km3/s2],[km] & [km]
end

for i = 1:size(pars.INPUTS.idMoon,2)
    pars.Moon.Vel(i)    = sqrt(pars.Planet.mu/pars.Moon.OrbRad(i));           %[km/s] Moon Orbital velocity
    pars.Moon.Period(i) = 2*pi*sqrt(pars.Moon.OrbRad(i)^3/pars.Planet.mu);    %[s] Moon orbital period
    pars.Moon.HillSph(i) = pars.Moon.OrbRad(i)*( pars.Moon.mu(i)/(3*(pars.Moon.mu(i) + pars.Planet.mu)))^(1/3);     %[km] Moon Hill's Sphere
end

% --- Parameters required for Linked Conic Calculation ---
pars.GroundTr.npoints      = 30e3; % N° of time steps for ground-track propagation
pars.GroundTr.t_prop       = 5;    %[minutes] Time of flyby hyperbola propagation 
pars.INPUTS.Flyby.hMapping = 300;  %[km] Max altitude to consider mapping 


%% ========================================================================
%  2. CALCULATE INITIAL STATES FOR PROPAGATION
%  ========================================================================

% Get the state of Enceladus at the time of the flyby
[r_enceladus_at_flyby, v_enceladus_at_flyby] = EphSS_car_spice2(602, pars.INPUTS.epoch0 , true, spiceParam);

% --- State for BACKWARD propagation (from nodein) ---
[vvinfin, r0_sc_in, v0_sc_in, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodein(1), nodein(2), nodein(3), pars.INPUTS.epoch0, pars.INPUTS.idMoon, pars.INPUTS.idCentral);

% --- State for FORWARD propagation (from nodeout) ---
[vvinfout, r0_sc_out, v0_sc_out, vvga] = vinfAlphaCrank_to_VinfCARTClaudio(nodeout(1), nodeout(2), nodeout(3), pars.INPUTS.epoch0, pars.INPUTS.idMoon, pars.INPUTS.idCentral );


% Determine maximum bending due to flyby
rp_flyby  = pars.INPUTS.Flyby.min_h + pars.Moon.EquRad;           %[km]
e_fly     = 1 + ((rp_flyby*pars.INPUTS.V_inf^2)/pars.Moon.mu);    %[-] 
delta_max = 2*asin(1/e_fly);                                      %[rad]
pars.delta_max = delta_max;

% check
[Flyby] = Flyby_BuildUp_claudio(nodein, nodeout, pars);

% coherence check 
stin = [r0_sc_in, v0_sc_in];
stout = [r0_sc_out, v0_sc_out];
errin = stin-Flyby.State_In;
errout = stout-Flyby.State_Out;
% r0_sc_in = offsetR0_sc(r0_sc_in, pars);
% r0_sc_out = offsetR0_sc(r0_sc_out, pars);

% --> convert in body-fixed reference frame (this is GTOC6 ref. frame)
b1        = -r0_sc_out./norm(r0_sc_out);
b3        = cross(r0_sc_out, vvga)./norm(cross(r0_sc_out, vvga));
b2        = cross(b3,b1);

[vvinfouBM, delta] = infVelBeforeDefect(vvinfin, vvinfout, delta_max, vvga);

Rm        = [ b1' b2' b3' ]';
vvinfin  =  [Rm*vvinfin']';
vvinfouBM = [Rm*vvinfouBM']';

% calculate pericentre state
Energy     = 0.5*norm(vvinfin)^2;
sma        = -mu_enceladus/(2*Energy);
ecc        = 1/(sin(delta/2));
rp         = sma*(1 - ecc);
hhat       = cross( vvinfin, vvinfouBM )./norm(cross( vvinfin, vvinfouBM ));
vp         = sqrt( norm(vvinfin)^2 + 2*mu_enceladus/rp );

% pericentre state, enceladus body fixed
rrp        = rp.*( vvinfin - vvinfouBM )./norm( vvinfin - vvinfouBM );
vvp        = vp.*cross( hhat, rrp./rp );

Rm_inv = Rm';

% Convert the pericenter position and velocity vectors back to the
% Saturn-centric reference frame.
rrp_saturn_centric = (Rm_inv * rrp')'  + r_enceladus_at_flyby;
vvp_saturn_centric = (Rm_inv * vvp')'  + v_enceladus_at_flyby;


% Ensure vectors are columns for the propagator
% r0_sc_in = r0_sc_in(:);
% v0_sc_in = v0_sc_in(:);
% r0_sc_out = r0_sc_out(:);
% v0_sc_out = v0_sc_out(:);

rrp_saturn_centric = rrp_saturn_centric(:);
vvp_saturn_centric = vvp_saturn_centric(:);

%% ========================================================================
%  3. PROPAGATE THE TRAJECTORIES
%  ========================================================================

% --- Setup Propagation Duration and Time Vectors ---
propagation_duration_days = 0.1;
duration_sec = propagation_duration_days * 86400;
time_steps = 30e3;

% Time vector for FORWARD propagation
time_vector_fwd = linspace(0, duration_sec, time_steps)';

% Time vector for BACKWARD propagation (negative duration)
time_vector_bwd = linspace(0, -duration_sec, time_steps)';

% --- Setup Ephemeris Handle ---
ephem_handle = @(t_sec) get_body_positions_wrapper(t_sec, pars.INPUTS.epoch0 , ...
    perturbingBodyNaifIDs, spiceParam);

% --- Run Forward Propagation from nodeout ---
[time_out_fwd, state_out_fwd] = propagateNBodyODE2(rrp_saturn_centric, vvp_saturn_centric, time_vector_fwd, ...
    mu_central_body, mu_TBs, ephem_handle, []);

% --- Run Backward Propagation from nodein ---
[time_out_bwd, state_out_bwd] = propagateNBodyODE2(rrp_saturn_centric, vvp_saturn_centric, time_vector_bwd, ...
    mu_central_body, mu_TBs, ephem_handle, []);

% --- Combine Trajectories ---
% backward propagation needs to be flipped to be in chronological order
time_out_bwd = flipud(time_out_bwd);
state_out_bwd = flipud(state_out_bwd);

% Combine the two trajectories
% Note: The pericenter point (t=0) will be duplicated, so remove one.
full_time_out = [time_out_bwd; time_out_fwd(2:end)];
full_state_out = [state_out_bwd; state_out_fwd(2:end, :)];
num_steps = length(full_time_out);

%% ========================================================================
%  4. ANIMATE 3D TRAJECTORY
%  ========================================================================
disp('Preparing data for animation...');

% --- First, calculate the full position history of Enceladus ---
r_enc_history = zeros(num_steps, 3);
for i = 1:num_steps
    current_mjd = pars.INPUTS.epoch0 + full_time_out(i) / 86400;
    [r_enc, ~] = EphSS_car_spice2(602, current_mjd, true, spiceParam);
    r_enc_history(i, :) = r_enc;
end

% --- Set up the animation figure ---
figure('Name', 'Animated Trajectory', 'Color', 'w');
ax = axes('Parent', gcf);
hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');
view(3);

xlabel(ax, 'ICRF X [km]');
ylabel(ax, 'ICRF Y [km]');
zlabel(ax, 'ICRF Z [km]');
title(ax, 'Spacecraft and Enceladus Trajectory');

% --- Initialize plot objects for animation for better performance ---
plot3(ax, r_enc_history(:,1), r_enc_history(:,2), r_enc_history(:,3), ...
    '--', 'Color', [0.5 0.5 0.5], 'DisplayName', 'Enceladus Full Path');
sc_path_plot = plot3(ax, NaN, NaN, NaN, 'b-', 'LineWidth', 2, 'DisplayName', 'Spacecraft');
enc_path_plot = plot3(ax, NaN, NaN, NaN, 'g-', 'LineWidth', 2, 'DisplayName', 'Enceladus');
sc_head_plot = plot3(ax, NaN, NaN, NaN, 'bo', 'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
enc_head_plot = plot3(ax, NaN, NaN, NaN, 'go', 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
legend(ax, 'show', 'Location', 'northeast');

% --- Animation Loop ---
disp('Starting animation...');
animation_step = 250; 
for i = 1:animation_step:num_steps
    set(sc_path_plot, 'XData', full_state_out(1:i, 1), 'YData', full_state_out(1:i, 2), 'ZData', full_state_out(1:i, 3));
    set(sc_head_plot, 'XData', full_state_out(i, 1), 'YData', full_state_out(i, 2), 'ZData', full_state_out(i, 3));
    set(enc_path_plot, 'XData', r_enc_history(1:i, 1), 'YData', r_enc_history(1:i, 2), 'ZData', r_enc_history(1:i, 3));
    set(enc_head_plot, 'XData', r_enc_history(i, 1), 'YData', r_enc_history(i, 2), 'ZData', r_enc_history(i, 3));
    drawnow;
end
disp('Animation finished.');

set(sc_path_plot, 'XData', full_state_out(:, 1), 'YData', full_state_out(:, 2), 'ZData', full_state_out(:, 3));
set(enc_path_plot, 'XData', r_enc_history(:, 1), 'YData', r_enc_history(:, 2), 'ZData', r_enc_history(:, 3));
set(sc_head_plot, 'XData', full_state_out(end, 1), 'YData', full_state_out(end, 2), 'ZData', full_state_out(end, 3));
set(enc_head_plot, 'XData', r_enc_history(end, 1), 'YData', r_enc_history(end, 2), 'ZData', r_enc_history(end, 3));
drawnow;


%% ========================================================================
%  5. COMPUTE N-BODY GROUND TRACK ON ENCELADUS 
%  ========================================================================
disp('Computing ground track from N-body propagation...');
lats = zeros(num_steps, 1);
longs = zeros(num_steps, 1);
alts = zeros(num_steps, 1);

for i = 1:num_steps
    current_mjd = pars.INPUTS.epoch0  + full_time_out(i) / 86400;
    sc_pos_icrf = full_state_out(i, 1:3); 
    r_enc = r_enc_history(i, :);
    r_iau_unit_vec = icrf2enceladus(sc_pos_icrf, r_enc, current_mjd);
    
    lon_rad_temp = atan2(r_iau_unit_vec(2), r_iau_unit_vec(1));
    longs(i) = lon_rad_temp + (lon_rad_temp < 0) * 2 * pi;
    lats(i) = asin(r_iau_unit_vec(3));
    
    sc_pos_relative = sc_pos_icrf - r_enc;
    alts(i) = norm(sc_pos_relative) - pars.Moon.EquRad; % Use radius from pars
end

lats_deg = rad2deg(lats);
longs_deg = rad2deg(longs); 

[~, idx_periapsis] = min(alts);
rp_lat_deg = lats_deg(idx_periapsis);
rp_long_deg = longs_deg(idx_periapsis);


%% ========================================================================
%  5.B COMPUTE LINKED CONIC REFERENCE GROUNDTRACK
%  ========================================================================
disp('Computing reference linked conic groundtrack...');

% NOTE: This section assumes you have the 'Flyby_BuildUp.m' function from 
% your demo script in your MATLAB path. This function calculates the 
% groundtrack using a linked-conic approximation.
[Flyby_LinkedConic] = Flyby_BuildUp(nodein, nodeout, pars);


%% ========================================================================
%  6. FILTER AND PLOT MAPPING GROUNDTRACK COMPARISON
%  ========================================================================

% --- Define the maximum mapping altitude  ---
mapping_altitude_km = pars.INPUTS.Flyby.hMapping;

% --- Find indices for the N-body track below the mapping threshold ---
valid_indices = alts <= mapping_altitude_km;

% --- Extract the relevant data for N-body mapping track ---
lats_for_mapping_rad = lats(valid_indices);
longs_for_mapping_rad = longs(valid_indices);

% --- Create a structure with the filtered data for the plotting function ---
NBody_Mapping_Data.lats = lats_for_mapping_rad;
NBody_Mapping_Data.longs = longs_for_mapping_rad;
NBody_Mapping_Data.rp_lat = rp_lat_deg;   % Periapsis in degrees
NBody_Mapping_Data.rp_long = rp_long_deg; % Periapsis in degrees

% --- Plot the results ---
figure('Name', 'Enceladus Ground Track Comparison', 'Color', [1 1 1]);
hold on;
plotTextureLatLong(pars.INPUTS.idMoon, pars.INPUTS.idCentral, 1); 
axis normal;
grid on;

% --- Plot N-Body Ground Track (Orange) ---
if isempty(lats_for_mapping_rad)
    disp(['No points on the N-body ground track were below mapping altitude of ', num2str(mapping_altitude_km), ' km.']);
else
    disp(['Plotting N-body track with ' num2str(numel(lats_for_mapping_rad)) ' points below ' num2str(mapping_altitude_km) ' km.']);
    Plot_Flyby_GT(NBody_Mapping_Data, [0.8500, 0.3250, 0.0980]); 
end

% --- Plot Linked Conic Ground Track (Blue) ---
% Note: The demo script used a loop for plotting. If Flyby_BuildUp returns
% a single structure for a single flyby, this direct call is sufficient.
if ~isempty(Flyby_LinkedConic)
     disp('Plotting reference linked conic ground track.');
     Plot_Flyby_GT(Flyby_LinkedConic, [0, 0.4470, 0.7410]); % Blue color
else
    disp('No linked conic ground track data was computed.');
end

% --- Add Title and Legend ---
title(['Ground Track Comparison (Altitude <= ' num2str(mapping_altitude_km) ' km)']);
% Create dummy plot handles for the legend
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'-', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2);
h(2) = plot(NaN,NaN,'-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2);
legend(h, 'N-Body Propagation', 'Linked Conic Approx.', 'Location', 'best');


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