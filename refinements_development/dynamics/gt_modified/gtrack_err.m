clc;
clear all;
close all;

%% 1. Define Inputs
% ------------------
% Flyby parameters
pars.INPUTS.idCentral = 6;      % Saturn
pars.INPUTS.idMoon    = 1;      % Enceladus flyby (NOTE: NAIF ID for Enceladus is 602. The code seems to use an internal index `1` for the moon)
pars.INPUTS.V_inf     = 4.0;    % km/s
pars.INPUTS.alpha     = 0.15;      % rad
pars.INPUTS.k         = 0.1;

% gtrack pars
pars.INPUTS.Flyby.min_h    = 25;   %[km] Minimum flyby altitude
pars.GroundTr.npoints      = 30e3; % N° of time steps for ground-track propagation
pars.GroundTr.t_prop       = 5;    %[minutes] Time of flyby hyperbola propagation 
pars.INPUTS.Flyby.hMapping = 300;  %[km] Max altitude to

% SPICE and Time parameters
% Make sure these kernel files are in your MATLAB path
kernels = {'sat441.bsp', 'naif0012.tls', 'pck00010.tpc', 'de440.bsp'};
loadSpiceKernels(kernels)
spiceParam.frame    = 'J2000';
spiceParam.abcorr   = 'NONE';
spiceParam.observer = '699'; % Saturn NAIF ID
t0_mjd = date2mjd2000([2025, 1, 1, 0, 0, 0]);

% N-Body Perturbers
% perturbingBodyNaifIDs = [602, 603, 605, 604, 601, 606, 10]; % Enceladus, Tethys, Rhea, Dione, Mimas, Titan, Sun
perturbingBodyNaifIDs = [606]; % Titan only

%% 2. Call the Function
[rp_target, rp_nbody] = calculatePericenterComparison(pars, spiceParam, t0_mjd, perturbingBodyNaifIDs, true);


%% 3. Display and Compare Results

% --- Display Target Pericenter ---
fprintf('--- Target Pericenter (Linked Conics) ---\n');
fprintf('Altitude (km):      %12.4f\n', rp_target.alt_km);
fprintf('Latitude (deg):     %12.4f\n', rp_target.lat_deg);
fprintf('Longitude (deg):    %12.4f\n\n', rp_target.lon_deg);

% --- Display N-Body Pericenter ---
fprintf('--- Actual Pericenter (N-Body Sim) ---\n');
fprintf('Altitude (km):      %12.4f\n', rp_nbody.alt_km);
fprintf('Latitude (deg):     %12.4f\n', rp_nbody.lat_deg);
fprintf('Longitude (deg):    %12.4f\n\n', rp_nbody.lon_deg);

% --- Calculate and Display Errors ---
alt_error_km  = rp_nbody.alt_km - rp_target.alt_km;
lat_error_deg = rp_nbody.lat_deg - rp_target.lat_deg;
lon_error_deg = rp_nbody.lon_deg - rp_target.lon_deg;

% Handle longitude wrap-around for a more accurate difference
if abs(lon_error_deg) > 180
    lon_error_deg = lon_error_deg - sign(lon_error_deg) * 360;
end

fprintf('--- Error (N-Body minus Target) ---\n');
fprintf('Altitude Error (km):  %+12.4f\n', alt_error_km);
fprintf('Latitude Error (deg): %+12.4f\n', lat_error_deg);
fprintf('Longitude Error (deg):%+12.4f\n', lon_error_deg);
fprintf('====================================================\n');

% Unload kernels at the end of the script
cspice_kclear;