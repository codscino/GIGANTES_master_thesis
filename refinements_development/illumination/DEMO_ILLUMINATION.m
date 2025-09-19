%% MAIN SCRIPT TO PLOT FLYBY GROUNDTRACKS WITH ILLUMINATION-BASED COLORING %%
clear all; close all; clc;

%% DEFINE PLANET, MOON PARAMETERS & CONSTANTS (COMMON TO BOTH SCENARIOS) %%
% Define Central Body & Moon of interest
pars.INPUTS.idCentral = 6; % (6 = Saturn)
pars.INPUTS.idMoon = [1]; % (1 = Enceladus)

% Define NAIF IDs for SPICE calls
pars.INPUTS.NAIFCentral = 699; % Saturn
pars.INPUTS.NAIFMoon = 602; % Enceladus

% Retrieve planetary constants
[pars.Planet.mu, pars.Planet.EquRad, pars.Planet.OrbRad, pars.Planet.hmin] = planetConstants(pars.INPUTS.idCentral); %[km3/s2],[km],[km] & [km]

% Retrieve Moon parameters
if pars.INPUTS.idCentral == 6
    [pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = satMoonsConstants(pars.INPUTS.idMoon); %[km],[km3/s2],[km] & [km]
else
    error('This script is configured for Saturn and its moons.');
end

% Calculate additional moon orbital properties
for i = 1:size(pars.INPUTS.idMoon,2)
    pars.Moon.Vel(i) = sqrt(pars.Planet.mu/pars.Moon.OrbRad(i)); %[km/s] Moon Orbital velocity
    pars.Moon.Period(i) = 2*pi*sqrt(pars.Moon.OrbRad(i)^3/pars.Planet.mu); %[s] Moon orbital period
    pars.Moon.HillSph(i) = pars.Moon.OrbRad(i)*( pars.Moon.mu(i)/(3*(pars.Moon.mu(i) + pars.Planet.mu)))^(1/3); %[km] Moon Hill's Sphere
end

% --- SPICE and Time Parameters ---
% Ensure you have these kernels in your MATLAB path or provide the full path
kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels);

spiceParam.frame = 'J2000';
spiceParam.abcorr = 'NONE';
spiceParam.observer = '699'; % Saturn Center NAIF ID

%% DEFINE FLYBY PARAMETERS & CONSTANTS (COMMON TO BOTH SCENARIOS) %%
% Define Flyby Hyperbolic Excess Velocity
pars.INPUTS.V_inf = 4; %[km/s]
vinfin = pars.INPUTS.V_inf;

% Define incoming/outgoing vector angles
alfain = 1;
kin = 0;
kou = 0.09;

% Define parameters regarding the flyby
pars.INPUTS.Flyby.min_h = 25; %[km] Minimum flyby altitude
pars.GroundTr.npoints = 3000; % N° of time steps for ground-track propagation
pars.GroundTr.t_prop = 5; %[minutes] Time of ONE-WAY hyperbola propagation (from pericenter)
pars.INPUTS.Flyby.hMapping = 300; %[km] Max altitude to consider for mapping

% Calculate the total duration of the flyby
% The total propagation time is twice the one-way time (before and after pericenter)
flyby_duration_minutes = pars.GroundTr.t_prop * 2; % [minutes]

% Determine maximum bending angle due to flyby
rp_flyby = pars.INPUTS.Flyby.min_h + pars.Moon.EquRad; %[km]
pars.INPUTS.Flyby.rp_flyby = rp_flyby;
e_fly = 1 + ((rp_flyby*pars.INPUTS.V_inf^2)/pars.Moon.mu); %[-]
pars.INPUTS.Flyby.e_fly = e_fly;
delta_max = 2*asin(1/e_fly); %[rad]
pars.delta_max = delta_max;

% Define incoming & outgoing nodes based on velocity and angles
nodein = [vinfin, alfain, kin]; %[km/s, rad, rad]
nodeout = [vinfin, alfain, kou];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- First Plot with eclipse colouring --- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('--- Processing First Request (2025) ---\n');
tic;

% --- Setup Time ---
date0_scen1 = [2025, 5, 10, 21, 46, 59];
% date0_scen1 = [2025, 5, 6, 19, 4, 34];
T0_scen1 = date2mjd2000(date0_scen1); % Convert to mjd2000
pars.INPUTS.epoch0 = T0_scen1;

% --- Compute Flyby Characteristics for the first scenario ---
[Flyby1] = Flyby_BuildUp_claudio(nodein, nodeout, pars);

% --- Setup Figure ---
fig1 = figure('Name', 'Flyby Groundtrack (2025)', 'Color', [1 1 1]);
hold on;

% Add date/time to title
dateVec1 = mjd20002date(T0_scen1);
dateStr1 = sprintf('%04d-%02d-%02d %02d:%02d:%02.0f', ...
    dateVec1(1), dateVec1(2), dateVec1(3), ...
    dateVec1(4), dateVec1(5), dateVec1(6));
title(sprintf('Enceladus Flyby Groundtrack with Illumination\nPericentre: %s', dateStr1));
xlabel('Longitude (deg)'); ylabel('Latitude (deg)');

% --- Plotting Parameters ---
colors1 = cool(length(Flyby1));
use_parallel = true;
coloured_illumination = true;
background_illumination = true;
step_deg = 2;

% --- Loop and Plot Each Flyby ---
for i = 1:size(Flyby1, 2)
    fprintf('Plotting Flyby Trajectory #%d with illumination...\n', i);
    
    Plot_Flyby_GT_illuminated(Flyby1(i), colors1(i,:), T0_scen1, pars,...
        flyby_duration_minutes, use_parallel, coloured_illumination,...
        background_illumination, step_deg);
end

xlim([0 360]); ylim([-90 90]);
hold off;
elapsed_time1 = toc;
fprintf('First plot completed in %.2f seconds\n\n', elapsed_time1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- Second plot with Incidence Angle Ranges --- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('--- Processing Second Request (2019) ---\n');
tic;

% --- Setup Time ---
date0_scen2 = [2019, 5, 10, 9, 37, 15];
T0_scen2 = date2mjd2000(date0_scen2); % Convert to mjd2000
pars.INPUTS.epoch0 = T0_scen2;

% --- Compute Flyby Characteristics for the second scenario ---
[Flyby2] = Flyby_BuildUp_claudio(nodein, nodeout, pars);

% --- Setup Figure ---
fig2 = figure('Name', 'Flyby Groundtrack with Incidence Angle (2019)', 'Color', [1 1 1]);
hold on;
plotTextureLatLong(pars.INPUTS.idMoon, pars.INPUTS.idCentral, 1);
axis normal; grid on;

% Add date/time to title
dateVec2 = mjd20002date(T0_scen2);
dateStr2 = sprintf('%04d-%02d-%02d %02d:%02d:%02.0f', ...
    dateVec2(1), dateVec2(2), dateVec2(3), ...
    dateVec2(4), dateVec2(5), dateVec2(6));
title(sprintf('Enceladus Flyby with Incidence Angle Constraints\nPericentre: %s', dateStr2));
xlabel('Longitude (deg)'); ylabel('Latitude (deg)');

% --- Plotting Parameters ---
colors2 = cool(length(Flyby2));
use_parallel = true;
coloured_illumination = true;
background_illumination = true;
step_deg = 0.5;
incidence_angle_ranges = [62,66; 68,90];

% --- Loop and Plot Each Flyby ---
for i = 1:size(Flyby2, 2)
    fprintf('Plotting Flyby Trajectory #%d with illumination...\n', i);

    Plot_Flyby_GT_illuminated(Flyby2(i), colors2(i,:), T0_scen2, pars,...
        flyby_duration_minutes, use_parallel, coloured_illumination,...
        background_illumination, step_deg, incidence_angle_ranges);
end

xlim([0 360]); ylim([-90 90]);
hold off;
elapsed_time2 = toc;
fprintf('Second plot completed in %.2f seconds\n', elapsed_time2);



