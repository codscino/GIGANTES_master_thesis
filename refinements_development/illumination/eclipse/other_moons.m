% Complete script to calculate the umbra and penumbra ratios of eclipses cast
% by Saturn and its moons on Enceladus, including eclipse duration and frequency.
% The umb_norm is the h_umb normalized with the Enceladus radius: 252.100 km
%
% This script includes:
% Part 1&2: Shadow size calculations
% Part 3: Eclipse duration calculations
% Part 4: Eclipse frequency (eclipse season) calculations

clear
clc

r_sun = 696340; % Sun radius in km

% Modified to capture mu_saturn
[mu_saturn, r_saturn, a_saturn, ~, ~, ~] = planetConstants(6);

%% Saturn umbra
[d_sat_enceladus, ~, r_enceladus, ~] = satMoonsConstants(1);
% Calculate Enceladus orbital period using Kepler's Third Law
T_enceladus = 2 * pi * sqrt(d_sat_enceladus^3 / mu_saturn); % Period in seconds

[h_umb_saturn, h_pen_saturn] = umbr_pen(r_sun, r_saturn, a_saturn, d_sat_enceladus);
ratio_sat = h_umb_saturn / r_enceladus;

% Print Enceladus radius
fprintf('Enceladus radius: %.3f km\n', r_enceladus);

% Calculate umb/pen ratio for Saturn
umb_pen_ratio_saturn = h_umb_saturn / h_pen_saturn;

% Now modify Saturn's shadow stats print line
fprintf('%-10s: h_umb = %7.3f km, umb_norm = %6.3f, umb/pen ratio = %6.3f\n', ...
        'Saturn', h_umb_saturn, ratio_sat, umb_pen_ratio_saturn);

%% for-loop moons (Parts 1&2: Shadow calculations)
% d_sm → saturn-moon distance
% r_m  → radius of the moon

% Define the moon names and IDs
moonNames = {'Mimas','Enceladus','Tethys','Dione','Rhea','Titan'};
moonIDs   = 0:5;            % [0=Mimas, 1=Enceladus, …, 5=Titan]

% Remove Enceladus (ID == 1) in one line:
mask      = moonIDs~=1;
moonNames = moonNames(mask);
moonIDs   = moonIDs(mask);

nMoons    = numel(moonIDs);
h_umb     = zeros(1,nMoons);
h_pen     = zeros(1,nMoons);
h_umb_norm= zeros(1,nMoons);
umb_pen_ratio = zeros(1,nMoons);
d_sm_all = zeros(1, nMoons); % Store semi-major axes for later use
r_m_all = zeros(1, nMoons); % Store radii for later use
T_moon_all = zeros(1, nMoons); % Store periods for later use

fprintf('\n--- Eclipse Potential by Saturn''s Moons on Enceladus ---\n');
for k = 1:nMoons
    id = moonIDs(k);
    [d_sm, ~, r_m, ~] = satMoonsConstants(id);
    d_sm_all(k) = d_sm;
    r_m_all(k) = r_m;

    % Calculate orbital period using Kepler's Third Law
    T_moon_all(k) = 2 * pi * sqrt(d_sm^3 / mu_saturn); % Period in seconds

    [h_umb(k), h_pen(k)] = umbr_pen(r_sun, r_m, a_saturn, abs(d_sm - d_sat_enceladus));
    h_umb_norm(k)   = h_umb(k) / r_enceladus;    % normalized to Enceladus radius
    umb_pen_ratio(k)= h_umb(k) / h_pen(k);       % umb/pen ratio
end

% Display
for k = 1:nMoons
    fprintf('%-10s: h_umb = %7.3f km, umb_norm = %6.3f, umb/pen ratio = %6.3f\n', ...
            moonNames{k}, h_umb(k), h_umb_norm(k), umb_pen_ratio(k));
end

%% Part 3: Eclipse Duration Calculation
fprintf('\n--- Part 3: Eclipse Duration Calculation ---\n');
fprintf('Calculating maximum eclipse durations based on relative orbital velocities.\n\n');

% Calculate Enceladus's orbital velocity (km/s)
v_enceladus = (2 * pi * d_sat_enceladus) / T_enceladus;

fprintf('%-10s | %-20s | %-25s\n', 'Moon', 'Relative Velocity (km/s)', 'Max Eclipse Duration (min)');
fprintf('%-10s-+-%-20s-+-%-25s\n', repmat('-',1,10), repmat('-',1,20), repmat('-',1,25));

for k = 1:nMoons
    % Get moon data
    [d_sm, ~, ~, ~] = satMoonsConstants(moonIDs(k));

    % Calculate moon's orbital period using Kepler's Third Law
    T_moon = 2 * pi * sqrt(d_sm^3 / mu_saturn); % Period in seconds

    % Calculate moon's orbital velocity (km/s)
    v_moon = (2 * pi * d_sm) / T_moon;

    % Calculate relative velocity
    v_rel = abs(v_enceladus - v_moon);

    % Calculate maximum eclipse duration in minutes
    % Duration = (2 * umbra_radius) / relative_velocity
    t_duration_min = (2 * h_umb(k)) / v_rel / 60;


    fprintf('%-10s | %20.3f | %25.2f\n', moonNames{k}, v_rel, t_duration_min);

end
%% Part 4: Eclipse Frequency Calculation
fprintf('\n--- Part 4: Eclipse Frequency Calculation ---\n');
fprintf('Calculating eclipse season windows around Saturn''s equinoxes.\n');
fprintf('(Eclipses are only possible during these windows due to orbital geometry)\n\n');

% Saturn's orbital and axial parameters
axial_tilt_deg = 26.73; % Saturn's axial tilt in degrees
T_saturn_orbit_days = 29.5 * 365.25; % Saturn's orbital period in Earth days

fprintf('Saturn''s axial tilt: %.2f degrees\n', axial_tilt_deg);
fprintf('Saturn''s orbital period: %.1f Earth years (%.0f days)\n', ...
        T_saturn_orbit_days/365.25, T_saturn_orbit_days);
fprintf('Equinoxes occur every ~%.1f years\n\n', T_saturn_orbit_days/365.25/2);

% Calculate the ring plane change rate at equinox
orbital_rate = 360 / T_saturn_orbit_days; % Saturn's mean motion (deg/day)
ring_plane_rate = orbital_rate * sind(axial_tilt_deg); % deg/day at equinox

fprintf('Ring plane angle change rate at equinox: %.5f deg/day\n', ring_plane_rate);
fprintf('This means the ring plane tilts by ~1 degree every %.1f days\n\n', 1/ring_plane_rate);

fprintf('%-10s | %-20s | %-30s\n', 'Moon', 'Critical Angle (deg)', 'Eclipse Season Window (days)');
fprintf('%-10s-+-%-20s-+-%-30s\n', repmat('-',1,10), repmat('-',1,20), repmat('-',1,30));

for k = 1:nMoons
    % Maximum allowable shadow displacement for eclipse to occur
    z_max = r_enceladus + h_umb(k);

    % Distance between the eclipsing moon and Enceladus
    d_moon_enceladus = abs(d_sm_all(k) - d_sat_enceladus);

    % Calculate critical angle (degrees)
    theta_crit_deg = atand(z_max / d_moon_enceladus);

    % Calculate eclipse season window (Earth days)
    % Window = time for ring plane to sweep through ±θ_crit
    season_window_days = 2 * theta_crit_deg / ring_plane_rate;

    % Display results
    fprintf('%-10s | %20.3f | %30.1f\n', ...
            moonNames{k}, theta_crit_deg, season_window_days);
end

fprintf('\nNote: Eclipse seasons occur twice per Saturn orbit (every ~%.1f years),\n', ...
        T_saturn_orbit_days/365.25/2);
fprintf('      centered around Saturn''s equinoxes.\n');
fprintf('      These are the windows when eclipses are geometrically possible.\n');

%% Part 5: Synodic Periods and Eclipse Opportunities
fprintf('\n--- Part 5: Synodic Periods and Eclipse Opportunities ---\n');
fprintf('Calculating how often moons align with Enceladus for potential eclipses.\n\n');

% Calculate the ring plane change rate for use in opportunity calculations
orbital_rate = 360 / T_saturn_orbit_days; % Saturn's mean motion (deg/day)
ring_plane_rate = orbital_rate * sind(axial_tilt_deg); % deg/day at equinox

fprintf('Ring plane angle change rate at equinox: %.5f deg/day\n', ring_plane_rate);
fprintf('This means the ring plane tilts by ~1 degree every %.1f days\n\n', 1/ring_plane_rate);

% First, display the orbital periods for reference
fprintf('Orbital Periods:\n');
fprintf('Enceladus: %.3f days\n', T_enceladus / 86400);
for k = 1:nMoons
    fprintf('%-10s: %.3f days\n', moonNames{k}, T_moon_all(k) / 86400);
end
fprintf('\n');

% Calculate synodic periods
fprintf('%-10s | %-20s | %-25s | %-30s\n', ...
        'Moon', 'Synodic Period (days)', 'Eclipse Window (days)', 'Potential Opportunities');
fprintf('%-10s-+-%-20s-+-%-25s-+-%-30s\n', ...
        repmat('-',1,10), repmat('-',1,20), repmat('-',1,25), repmat('-',1,30));

for k = 1:nMoons
    % Calculate synodic period with Enceladus
    % T_synodic = (T1 * T2) / |T2 - T1|
    T_synodic_seconds = abs((T_enceladus * T_moon_all(k)) / (T_moon_all(k) - T_enceladus));
    T_synodic_days = T_synodic_seconds / 86400; % Convert to days

    % Recalculate the eclipse window with the correct formula
    z_max = r_enceladus + h_umb(k);
    d_moon_enceladus = abs(d_sm_all(k) - d_sat_enceladus);
    theta_crit_deg = atand(z_max / d_moon_enceladus);

    % Use the more accurate formula for small angles
    season_window_days = 2 * theta_crit_deg / ring_plane_rate;

    % Calculate number of potential eclipse opportunities(divided by 2
    % beacuse it only happens at conjunctions)
    % This is how many times the moons align during the eclipse window
    num_opportunities = (season_window_days / T_synodic_days) / 2;

    % Display results

        if num_opportunities < 1
            fprintf('%-10s | %20.3f | %25.2f | %30s\n', ...
                    moonNames{k}, T_synodic_days, season_window_days, '<1 (may miss)');
        else
            fprintf('%-10s | %20.3f | %25.2f | %30.1f\n', ...
                    moonNames{k}, T_synodic_days, season_window_days, num_opportunities);
        end
end