function [rrf, vvf, kep2] = FGKepler_dt_robust(kep1, dt, mu)
% DESCRIPTION
% Robust Kepler propagation with proper angle wrapping handling.
% This version prevents discontinuities in orbital motion.
%
% INPUT
% - kep1 : 1x6 vector with initial keplerian elements (see car2kep.m)
% - dt : time of flight [s]
% - mu : gravitational parameter of the central body [km3/s2]
%
% OUTPUT
% - rrf : 1x3 vector of spacecraft position after the propagation [km]
% - vvf : 1x3 vector of spacecraft velocity after the propagation [km/s]
% - kep2 : 1x6 vector of final keplerian elements after the propagation
%
% -------------------------------------------------------------------------

% Calculate mean motion
if kep1(2) < 1
    n = sqrt(mu/kep1(1)^3);
elseif kep1(2) > 1
    n = sqrt(mu/(-kep1(1)^3));
else
    error('Parabolic orbits (e=1) not supported');
end

% Convert initial true anomaly to mean anomaly
M1 = theta2M(kep1(6), kep1(2));

% Propagate mean anomaly (this naturally handles wrapping)
M2 = M1 + n*dt;

% Normalize mean anomaly to [0, 2π] range to prevent numerical issues
M2 = mod(M2, 2*pi);

% Convert back to true anomaly
th2 = M2theta(M2, kep1(2));

% Ensure true anomaly is in [0, 2π] range
th2 = mod(th2, 2*pi);

% Build final Keplerian elements
kep2 = [kep1(1) kep1(2) kep1(3) kep1(4) kep1(5) th2];

% Convert to Cartesian coordinates
car2 = kep2car(kep2, mu);
rrf = car2(1:3);
vvf = car2(4:6);

end

% -------------------------------------------------------------------------
% Alternative: Direct circular orbit propagation for better numerical stability
% -------------------------------------------------------------------------

function [rrf, vvf, kep2] = propagateCircularOrbit(kep1, dt, mu)
% For circular orbits (e ≈ 0), use direct trigonometric propagation
% This avoids conversion issues between mean and true anomaly

a = kep1(1);           % Semi-major axis
inc = kep1(3);         % Inclination
Omega = kep1(4);       % RAAN
omega = kep1(5);       % Argument of periapsis
theta0 = kep1(6);      % Initial true anomaly

% Mean motion for circular orbit
n = sqrt(mu/a^3);

% For circular orbits, mean anomaly = true anomaly
theta_final = theta0 + n*dt;

% Keep angle in reasonable range (unwrapped)
% Don't use mod() here to maintain continuity
kep2 = [a, 0, inc, Omega, omega, theta_final];

% Convert to Cartesian
car2 = kep2car(kep2, mu);
rrf = car2(1:3);
vvf = car2(4:6);

end

% -------------------------------------------------------------------------
% Modified Saturn moons ephemeris function
% -------------------------------------------------------------------------

function [rr, vv, kep] = approxEphemSatMoons_cc_robust(idmoon, t)
% Robust version of Saturn moons ephemeris with better angle handling

muidcentral = 37931005.114; % Saturn's gravitational parameter
tref = date2mjd2000([2030 1 1 0 0 0]); % Reference epoch

% Initial Keplerian elements at reference time
if idmoon == 1 % Enceladus
    kep0 = [237948, 0, 0, 0, 0, deg2rad(250.0815)]; % Simplified initial true anomaly
elseif idmoon == 2 % Tethys
    kep0 = [294619, 0, 0, 0, 0, deg2rad(1.9327)];
elseif idmoon == 3 % Dione
    kep0 = [377396, 0, 0, 0, 0, deg2rad(338.7030)];
elseif idmoon == 4 % Rhea
    kep0 = [527108, 0, 0, 0, 0, deg2rad(31.3014)];
elseif idmoon == 5 % Titan
    kep0 = [1221870, 0, 0, 0, 0, deg2rad(213.3514)];
end

dt = t - tref;

% Use circular orbit propagation for better stability
[rr, vv, kep] = propagateCircularOrbit(kep0, dt*86400, muidcentral);

end