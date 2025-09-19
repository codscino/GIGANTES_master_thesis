function [rr, vv, kep] = approxEphemSatMoons_cc_numerical(idmoon, t)
% DESCRIPTION
% Approximate ephemerides of Saturn moons using numerical integration
% instead of analytical Kepler propagation. This eliminates angle wrapping
% discontinuities that can cause apparent motion reversals.
%
% INPUT
% - idmoon : ID of the moon (see constants.m)
% - t : epoch at which the ephemerides are computed [MJD2000]
%
% OUTPUT
% - rr : 1x3 vector with moon position [km]
% - vv : 1x3 vector with moon velocity [km]
% - kep : 1x6 vector with keplerian elements (see car2kep.m)
%
% -------------------------------------------------------------------------

muidcentral = 37931005.114; % Saturn's gravitational parameter [km^3/s^2]
tref = date2mjd2000([2030 1 1 0 0 0]); % Reference epoch (MJD2000)

% Initial Keplerian elements at reference epoch
if idmoon == 1 % Enceladus
    kep0 = [237948, 0, 0, 0, 0, deg2rad(250.0815419418998)];
elseif idmoon == 2 % Tethys  
    kep0 = [294619, 0, 0, 0, 0, deg2rad(1.932732150115001)];
elseif idmoon == 3 % Dione
    kep0 = [377396, 0, 0, 0, 0, deg2rad(338.7029874381438)];
elseif idmoon == 4 % Rhea
    kep0 = [527108, 0, 0, 0, 0, deg2rad(31.30144905250815)];
elseif idmoon == 5 % Titan
    kep0 = [1221870, 0, 0, 0, 0, deg2rad(213.3513863591349)];
else
    error('Unknown moon ID: %d', idmoon);
end

% Time difference from reference epoch
dt = (t - tref) * 86400; % Convert from days to seconds

% Convert initial Keplerian elements to Cartesian coordinates
car0 = kep2car(kep0, muidcentral);
r0 = car0(1:3)';  % Initial position [3x1]
v0 = car0(4:6)';  % Initial velocity [3x1]

% Create time vector for integration
% For small time steps, we can integrate directly
% For efficiency, we'll use a single integration step
if abs(dt) < 1e-6
    % Very small time difference - just return initial state
    rr = r0';
    vv = v0';
    kep = kep0;
    return;
end

% Use propagateKeplerODE2 to integrate from t=0 to t=dt
time_vector = [0, dt]; % Start at 0, end at dt

% Integrate the orbit
[~, state_history] = propagateKeplerODE2(r0, v0, time_vector, muidcentral);

% Extract final state (last row of state_history)
final_state = state_history(end, :);
rr = final_state(1:3);  % Final position [1x3]
vv = final_state(4:6);  % Final velocity [1x3]

% Convert back to Keplerian elements for output
try
    kep = car2kep([rr, vv], muidcentral);
catch
    % In case of numerical issues, return the propagated Cartesian state
    % and set keplerian elements to NaN
    kep = [NaN, NaN, NaN, NaN, NaN, NaN];
    warning('Could not convert final state to Keplerian elements for moon %d', idmoon);
end

end

% -------------------------------------------------------------------------
% Alternative version with adaptive time stepping for better accuracy
% -------------------------------------------------------------------------

function [rr, vv, kep] = approxEphemSatMoons_cc_adaptive(idmoon, t)
% Version with adaptive time stepping for very long propagation times

muidcentral = 37931005.114; % Saturn's gravitational parameter
tref = date2mjd2000([2030 1 1 0 0 0]); % Reference epoch

% Initial Keplerian elements at reference epoch
if idmoon == 1 % Enceladus
    kep0 = [237948, 0, 0, 0, 0, deg2rad(250.0815419418998)];
    orbital_period = 2*pi*sqrt(kep0(1)^3/muidcentral); % seconds
elseif idmoon == 2 % Tethys
    kep0 = [294619, 0, 0, 0, 0, deg2rad(1.932732150115001)];
    orbital_period = 2*pi*sqrt(kep0(1)^3/muidcentral);
elseif idmoon == 3 % Dione
    kep0 = [377396, 0, 0, 0, 0, deg2rad(338.7029874381438)];
    orbital_period = 2*pi*sqrt(kep0(1)^3/muidcentral);
elseif idmoon == 4 % Rhea
    kep0 = [527108, 0, 0, 0, 0, deg2rad(31.30144905250815)];
    orbital_period = 2*pi*sqrt(kep0(1)^3/muidcentral);
elseif idmoon == 5 % Titan
    kep0 = [1221870, 0, 0, 0, 0, deg2rad(213.3513863591349)];
    orbital_period = 2*pi*sqrt(kep0(1)^3/muidcentral);
end

dt = (t - tref) * 86400; % Convert to seconds

% Convert to Cartesian
car0 = kep2car(kep0, muidcentral);
r0 = car0(1:3)';
v0 = car0(4:6)';

% Handle very small time differences
if abs(dt) < 1e-6
    rr = r0';
    vv = v0';
    kep = kep0;
    return;
end

% For long propagation times, use intermediate steps to maintain accuracy
max_step_size = orbital_period / 100; % 1% of orbital period per step
num_steps = max(2, ceil(abs(dt) / max_step_size));
time_vector = linspace(0, dt, num_steps);

% Propagate using ODE integration
[~, state_history] = propagateKeplerODE2(r0, v0, time_vector, muidcentral);

% Extract final state
final_state = state_history(end, :);
rr = final_state(1:3);
vv = final_state(4:6);

% Convert to Keplerian elements
try
    kep = car2kep([rr, vv], muidcentral);
catch
    kep = [NaN, NaN, NaN, NaN, NaN, NaN];
    warning('Keplerian conversion failed for moon %d at t=%.3f', idmoon, t);
end

end

% -------------------------------------------------------------------------
% Modified main ephemeris function to use numerical integration
% -------------------------------------------------------------------------

function [rr, vv, kep] = approxEphem_CC_numerical(idmoon, t, idcentral)
% DESCRIPTION
% Modified approximate ephemeris function using numerical integration
% for Saturn moons to avoid angle wrapping issues
%
% INPUT
% - idmoon : ID of the flyby body
% - t : epoch [MJD2000]  
% - idcentral : ID of the central body
%
% OUTPUT
% - rr : 1x3 position vector [km]
% - vv : 1x3 velocity vector [km/s]
% - kep : 1x6 keplerian elements

if idcentral == 1 % Solar System planets
    % Use original function (not implemented here)
    error('Solar System planets not implemented in numerical version');
elseif idcentral == 5 % Jupiter moons
    % Use original function (not implemented here)  
    error('Jupiter moons not implemented in numerical version');
elseif idcentral == 6 % Saturn moons
    [rr, vv, kep] = approxEphemSatMoons_cc_numerical(idmoon, t);
elseif idcentral == 7 % Uranus moons
    % Use original function (not implemented here)
    error('Uranus moons not implemented in numerical version');
end

end