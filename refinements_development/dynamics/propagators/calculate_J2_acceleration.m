function a_J2 = calculate_J2_acceleration(r_vec, mu, J2, r_eq)
% Computes the acceleration due to the J2 perturbation.
%
%
% Formula:
% a_J2 = -3/2 * J2 * (mu/r^2) * (r_eq/r)^2 * [ (1-5(z/r)^2)*x/r;
%                                             (1-5(z/r)^2)*y/r;
%                                             (3-5(z/r)^2)*z/r ]
%
% INPUTS:
%   r_vec   (3x1 double) The position vector of the spacecraft relative to
%           the central body [x; y; z] in km.
%   mu      (double) Gravitational parameter of the central body (km^3/s^2).
%   J2      (double) J2 zonal harmonic coefficient of the central body (dimensionless).
%   r_eq    (double) Equatorial radius of the central body (km).
%
% OUTPUTS:
%   a_J2    (3x1 double) The perturbing acceleration vector due to J2 (km/s^2).

if nargin < 2 % it is Saturn J2
    mu = getAstroConstants('Saturn', 'Mu');
    J2 = 16290.573/10^6; % source: https://atmos.nmsu.edu/data_and_services/atmospheres_data/Cassini/logs/Table%201Measured%20gravity%20harmonic%20coefficients.pdf
    r_eq = 60268; % equatorial Saturn radius(values from Spice)
end

% r_vec components
x = r_vec(1);
y = r_vec(2);
z = r_vec(3);

r_mag = norm(r_vec); %magnitude r_vec

% term appears frequently
z_over_r_sq = (z / r_mag)^2; 

% Calculate the main scalar coefficient of the formula
% This is the part: -3/2 * J2 * (mu/r^2) * (r_eq/r)^2
scalar_term = -1.5 * J2 * (mu / r_mag^2) * (r_eq / r_mag)^2;

% Calculate each component of the acceleration vector
ax = scalar_term * (1 - 5 * z_over_r_sq) * (x / r_mag);
ay = scalar_term * (1 - 5 * z_over_r_sq) * (y / r_mag);
az = scalar_term * (3 - 5 * z_over_r_sq) * (z / r_mag);

% Assemble the final acceleration vector
a_J2 = [ax; ay; az];

end