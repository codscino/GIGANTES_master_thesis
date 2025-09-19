%% Main Script to Analyze Enceladus Lighting Conditions
clc;
clear all;

%% Initialization
% Load SPICE kernels
kernels = {'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'}; 

loadSpiceKernels(kernels);

% Get Enceladus mean radius
radii_enc = cspice_bodvrd('602', 'RADII', 3);
R_enc = mean(radii_enc); % Enceladus average radius in km (252.1 km)


%% light dist
% T_mid = date2mjd2000([2025, 05, 10, 18, 0, 0]);
% 
% spiceParam.frame    = 'J2000';
% spiceParam.abcorr   = 'NONE';
% spiceParam.observer = '602';
% r_enc_mid = EphSS_car_spice(699, T_mid, true, kernels, spiceParam);
% 
% kernels = {'sat441.bsp','naif0012.tls'};
% 
% spiceParam1.frame    = 'J2000';
% spiceParam1.abcorr   = 'LT+S';
% spiceParam1.observer = '602';
% 
% r_enc_mid_lt = EphSS_car_spice(699, T_mid, true, kernels, spiceParam1);
% 
% aa = r_enc_mid_lt-r_enc_mid;

%% Create a Sample Data Table
T_start = date2mjd2000([2025, 05, 10, 21, 40, 0]);

% Define the duration and number of steps for the simulation
num_hours = 0.2;
num_steps = 3600; % A data point every 2 minutes

% Calculate the end time in MJD2000
duration_in_days = num_hours / 24;
T_end = T_start + duration_in_days;

% Generate the array of time steps using linspace for efficiency.
% The transpose (') makes it a column vector, matching the other arrays.
T_array = linspace(T_start, T_end, num_steps)';

% Define the fixed point on Enceladus's surface
fixed_lat_rad = deg2rad(0);
fixed_lon_rad = deg2rad(0);

% Create the corresponding latitude and longitude arrays
lat_array = ones(num_steps, 1) * fixed_lat_rad;
lon_array = ones(num_steps, 1) * fixed_lon_rad;

% Create the final MATLAB table
dataTable = table(lat_array, lon_array, T_array, 'VariableNames', {'lat', 'lon', 'T'});


%% Run the Analysis
% This single function call performs all the requested logic.
spiceParam.frame    = 'J2000';
spiceParam.abcorr   = 'LT+S';
spiceParam.observer = '602';

kernels = {'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'}; 

tic
lighting_status = calculate_lighting_conditions(dataTable, R_enc, spiceParam);
toc

% Append results to the table for easy viewing
dataTable.Lighting = lighting_status;


%% Post-Analysis: Find the First Penumbra Event and Compare to Reference

% The 'Lighting' column is a cell array of strings. We need to find the
% first row where the value is 'Penumbra'.
is_penumbra = strcmp(dataTable.Lighting, 'Penumbra');
first_penumbra_index = find(is_penumbra, 1, 'first');

% Check if a penumbra event was found and display the results.
if ~isempty(first_penumbra_index)
    % Extract the MJD2000 time from the dataTable at that index
    T_first_penumbra = dataTable.T(first_penumbra_index);
    
    % Convert the MJD2000 time back to a human-readable date vector
    date_vec_penumbra = mjd20002date(T_first_penumbra);
    
    % --- Validation against Stellarium ---
    % Define the reference time from Stellarium and convert it to MJD2000
    stellarium_date = date2mjd2000([2025, 05, 10, 21, 45, 23]);

    % from spice (at center of enceladus) 21:46, 55s
    % from spice (at surface point) 21:46, 51s
    % my value with UTC and 0.2h 21:48, 2s
    % my value with ET and 0.2h 21:46, 52s
    % my value with ET and 1h 21:46, 54s
    
    % Calculate the absolute time difference in seconds
    error_date_seconds = abs(stellarium_date - T_first_penumbra) * 86400;

    % --- Print the complete results ---
    fprintf('\n--- First Penumbra Event Analysis ---\n');
    fprintf('The first "Penumbra" event occurs at table index: %d\n', first_penumbra_index);
    fprintf('Calculated Event Time (UTC): %s\n', datestr(date_vec_penumbra, 'yyyy-mm-dd HH:MM:SS.FFF'));
    fprintf('Calculated Event Time (MJD2000): %.6f\n', T_first_penumbra);
    fprintf('Time Difference from Stellarium: %.3f seconds\n', error_date_seconds);
    
else
    % This block runs if no 'Penumbra' strings were found in the entire column
    fprintf('\nNo "Penumbra" events were found in the simulation time frame.\n');
end


