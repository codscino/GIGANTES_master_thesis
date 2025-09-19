% Setup code to generate data for Enceladus comparison animation
clc
clear all

% Load SPICE kernels
kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels);
spiceParam.abcorr = 'NONE';
spiceParam.observer = '699';  % Saturn
spiceParam.frame = 'J2000';

% Time setup - shorter time span for Enceladus orbit (about 1.37 days)
L = 100000;  % More points for smoother animation
date1 = date2mjd2000([2030 1 1 0 0 0]);
date2 = date2mjd2000([2041 1 1 0 0 0]);  % 2 days to see orbital motion
dates = linspace(date1, date2, L);

% Initialize arrays
r_enc_2body = zeros(L, 3);  % 2-body approximation positions
r_enc_nbody = zeros(L, 3);  % N-body SPICE positions

% Calculate Enceladus positions using both methods
fprintf('Calculating Enceladus positions...\n');
for i = 1:L
    if mod(i, 20) == 0
        fprintf('Progress: %d/%d\n', i, L);
    end
    
    % N-body SPICE calculation (Enceladus = 602)
    [rr_spice, ~, ~] = EphSS_car_spice2(602, dates(i), true, spiceParam);
    r_enc_nbody(i,:) = rr_spice;
    
    % 2-body approximation (assuming Enceladus corresponds to body 1 in your system)
    % You might need to adjust the body number based on your implementation
    [rr_2body, ~, ~] = approxEphem_CC(1, dates(i), 6);  % Adjust body number as needed
    r_enc_2body(i,:) = rr_2body;
end

% Create time vector in seconds (relative to start)
full_time_out_nb = (dates - dates(1)) * 86400; % Convert days to seconds

% Parameters structure (you'll need to define this based on your system)
pars.Moon.OrbRad = 238000;      % Enceladus orbital radius in km
pars.Moon.EquRad = 252;         % Enceladus radius in km
pars.EncPlotSize = 20;         % Size factor for plotting

% Call the animation function
liveplot_enceladus_comparison(r_enc_2body, r_enc_nbody, full_time_out_nb, dates, pars);

% Print summary statistics
diff_positions = r_enc_nbody - r_enc_2body;
diff_magnitudes = vecnorm(diff_positions', 2);

fprintf('\n=== Enceladus Position Comparison Statistics ===\n');
fprintf('Time span: %.2f days\n', (dates(end) - dates(1)));
fprintf('Number of data points: %d\n', L);
fprintf('Position differences:\n');
fprintf('  Mean difference: %.3f km\n', mean(diff_magnitudes));
fprintf('  Max difference: %.3f km\n', max(diff_magnitudes));
fprintf('  RMS difference: %.3f km\n', sqrt(mean(diff_magnitudes.^2)));
fprintf('  Min difference: %.3f km\n', min(diff_magnitudes));

% Optional: Save the data for later use
% save('enceladus_comparison_data.mat', 'r_enc_2body', 'r_enc_nbody', 'full_time_out_nb', 'dates', 'pars');