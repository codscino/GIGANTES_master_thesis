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


%% ========================================================================
%  2. CALCULATE INITIAL STATES FOR PROPAGATION
%  ========================================================================

% Get the state of Enceladus at the time of the flyby
[r_enceladus_at_flyby, v_enceladus_at_flyby] = EphSS_car_spice2(602, pars.INPUTS.epoch0 , true, spiceParam);

% --- State for BACKWARD propagation (from nodein) ---
[~, r0_sc_in, v0_sc_in, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodein(1), nodein(2), nodein(3), pars.INPUTS.epoch0, pars.INPUTS.idMoon, pars.INPUTS.idCentral);

% --- State for FORWARD propagation (from nodeout) ---
[~, r0_sc_out, v0_sc_out, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodeout(1), nodeout(2), nodeout(3), pars.INPUTS.epoch0, pars.INPUTS.idMoon, pars.INPUTS.idCentral );


% Determine maximum bending due to flyby
rp_flyby  = pars.INPUTS.Flyby.min_h + pars.Moon.EquRad;           %[km]
e_fly     = 1 + ((rp_flyby*pars.INPUTS.V_inf^2)/pars.Moon.mu);    %[-] 
delta_max = 2*asin(1/e_fly);                                      %[rad]
pars.delta_max = delta_max;

% Define parameters regarding the flyby
pars.INPUTS.Flyby.min_h    = 25;   %[km] Minimum flyby altitude
pars.GroundTr.npoints      = 30e3; % N° of time steps for ground-track propagation
pars.GroundTr.t_prop       = 5;    %[minutes] Time of flyby hyperbola propagation 
pars.INPUTS.Flyby.hMapping = 300;  %[km] Max altitude to consider mapping 

[Flyby] = Flyby_BuildUp_claudio(nodein, nodeout, pars);

% coherence check 
stin = [r0_sc_in, v0_sc_in];
stout = [r0_sc_out, v0_sc_out];

errin = stin-Flyby.State_In;
errout = stout-Flyby.State_Out;

r0_sc_in = offsetR0_sc(r0_sc_in, pars);
r0_sc_out = offsetR0_sc(r0_sc_out, pars);

% Ensure vectors are columns for the propagator
r0_sc_in = r0_sc_in(:);
v0_sc_in = v0_sc_in(:);
r0_sc_out = r0_sc_out(:);
v0_sc_out = v0_sc_out(:);

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
[time_out_fwd, state_out_fwd] = propagateNBodyODE2(r0_sc_out, v0_sc_out, time_vector_fwd, ...
    mu_central_body, mu_TBs, ephem_handle, []);

% --- Run Backward Propagation from nodein ---
[time_out_bwd, state_out_bwd] = propagateNBodyODE2(r0_sc_in, v0_sc_in, time_vector_bwd, ...
    mu_central_body, mu_TBs, ephem_handle, []);

% --- Combine Trajectories ---
% backward propagation needs to be flipped to be in chronological order
time_out_bwd = flipud(time_out_bwd);
state_out_bwd = flipud(state_out_bwd);

% Combine the two trajectories
% Note: The pericenter point (t=0) will be duplicated, so remove one.
full_time_out = [time_out_bwd; time_out_fwd(2:end)];
full_state_out = [state_out_bwd; state_out_fwd(2:end, :)];

%% ========================================================================
%  4. COMPUTE GROUND TRACK ON ENCELADUS
%  ========================================================================
num_steps = length(full_time_out);
lats = zeros(num_steps, 1);
longs = zeros(num_steps, 1);
alts = zeros(num_steps, 1);

for i = 1:num_steps
    current_mjd = pars.INPUTS.epoch0  + full_time_out(i) / 86400;
    
    sc_pos_icrf = full_state_out(i, 1:3)'; 
    [r_enc, ~] = EphSS_car_spice2(pars.INPUTS.idMoon, current_mjd, true, spiceParam);
    
    r_iau_unit_vec = icrf2enceladus(sc_pos_icrf, r_enc, current_mjd);
    
    % --- Calculate Latitude and Longitude from the body-fixed unit vector ---
    lon_rad_temp = atan2(r_iau_unit_vec(2), r_iau_unit_vec(1));   % range of [-pi, pi]
    
    % Convert negative angles to the [0, 2*pi] range
    if lon_rad_temp < 0
        longs(i) = lon_rad_temp + 2*pi;
    else
        longs(i) = lon_rad_temp;
    end
    
    % Latitude is straightforward 
    lats(i) = asin(r_iau_unit_vec(3));

    % --- Calculate Altitude separately ---
    sc_pos_relative = sc_pos_icrf - r_enc(:);
    alts(i) = norm(sc_pos_relative) - 252;
end

% Convert to degrees for plotting. Longitude will now correctly be in [0, 360].
lats_deg = rad2deg(lats);
longs_deg = rad2deg(longs); 

% Find the periapsis point to mark it
[~, idx_periapsis] = min(alts);
rp_lat_deg = lats_deg(idx_periapsis);
rp_long_deg = longs_deg(idx_periapsis);


%% ========================================================================
%  5. FILTER AND PLOT MAPPING GROUNDTRACK
%  ========================================================================

% --- Define the maximum mapping altitude  ---
mapping_altitude_km = pars.INPUTS.Flyby.hMapping;

% --- Find the indices where the altitude is below the mapping threshold ---
% This creates a logical array where '1' corresponds to a valid altitude.
valid_indices = alts <= mapping_altitude_km;

% --- Use the logical index to extract the relevant data for mapping ---
lats_for_mapping_rad = lats(valid_indices);
longs_for_mapping_rad = longs(valid_indices);

% --- Create a structure with the filtered data for the plotting function ---
Flyby_Mapping_Data.lats = lats_for_mapping_rad;
Flyby_Mapping_Data.longs = longs_for_mapping_rad;
Flyby_Mapping_Data.rp_lat = rp_lat_deg;   % Periapsis in degrees
Flyby_Mapping_Data.rp_long = rp_long_deg; % Periapsis in degrees

% --- Plot the results ---
figure('Name', 'Enceladus Mapping Ground Track', 'Color', [1 1 1]);
hold on;
plotTextureLatLong(pars.INPUTS.idMoon, pars.INPUTS.idCentral, 1); 
axis normal;
grid on;

% --- Check if any points are within the mapping altitude before plotting ---
if isempty(lats_for_mapping_rad)
    disp(['No points on the ground track were below the specified mapping altitude of ', num2str(mapping_altitude_km), ' km.']);
    title(['No Mapping Track Found (Altitude <= ' num2str(mapping_altitude_km) ' km)']);
else
    disp(['Plotting ' num2str(numel(lats_for_mapping_rad)) ' points below ' num2str(mapping_altitude_km) ' km.']);
    Plot_Flyby_GT(Flyby_Mapping_Data, [0.8500, 0.3250, 0.0980]); % Orange color
    title(['Ground Track for Mapping (Altitude <= ' num2str(mapping_altitude_km) ' km)']);
end


%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================
function positions = get_body_positions_wrapper(t_sec, t0_mjd_start, naif_ids, sp_param)
    current_mjd = t0_mjd_start + t_sec / 86400;
    num_bodies = length(naif_ids);
    positions = zeros(3, num_bodies);
    
    for k = 1:num_bodies
        % isSpice is set to true to use the SPICE-based part of your function
        [r, ~] = EphSS_car_spice2(naif_ids(k), current_mjd, true, sp_param);
        positions(:, k) = r(:); % Ensure column vector
    end
end

