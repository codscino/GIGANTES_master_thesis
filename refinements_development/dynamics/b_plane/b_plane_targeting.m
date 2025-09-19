function [B_R, B_T, B_vec] = b_plane_targeting(r_vec, v_vec, v_inf_sq, mu)
% Calculates B-plane parameters from a state vector.
%
% From algorithm 78 of Vallado
%
% Inputs:
%   r_vec       (3x1) Position vector of the spacecraft [km]
%   v_vec       (3x1) Velocity vector of the spacecraft [km/s]
%   v_inf_sq    (1x1) Square of the hyperbolic excess velocity (v_inf^2) [km^2/s^2]
%   mu          (1x1) Gravitational parameter of the central body [km^3/s^2]
%
% Outputs:
%   B_R         (1x1) Component of the B-vector along the R-axis [km]
%   B_T         (1x1) Component of the B-vector along the T-axis [km]
%   B_vec       (3x1) The B-vector itself [km] (optional output)
%
% Coordinate System Reference:
%   The algorithm assumes a standard inertial reference frame (e.g., ECI).
%   The K_hat vector [0; 0; 1] is the reference direction, typically the
%   pole of the reference plane (e.g., the celestial pole for the equator).
%

% --- Input Validation ---
% Ensure inputs are column vectors for consistent matrix operations
r_vec = r_vec(:);
v_vec = v_vec(:);

% --- Step 1: Calculate Angular Momentum and Eccentricity Vectors ---

% Specific angular momentum vector h
h_vec = cross(r_vec, v_vec);
h_mag = norm(h_vec);
h_hat = h_vec / h_mag;

% Eccentricity vector e
r_mag = norm(r_vec);
v_sq = dot(v_vec, v_vec);
e_vec = ((v_sq - mu/r_mag) * r_vec - dot(r_vec, v_vec) * v_vec) / mu;
e_mag = norm(e_vec);

% --- Step 2: Calculate Hyperbolic Orbit Parameters ---

% Semi-major axis a
a = -mu / v_inf_sq;

% Semi-minor axis b (magnitude of the B-vector)
% The negative sign is a convention used in this algorithm.
b = -a * sqrt(e_mag^2 - 1);

% Angle of the incoming asymptote phi_s
phi_s = acos(1 / e_mag);

% --- Step 3: Define the B-Plane Coordinate System (S, T, R) ---

% K_hat is the fundamental reference direction (e.g., Z-axis of ECI frame)
K_hat = [0; 0; 1];

% S_hat is the unit vector along the incoming asymptote
% S_hat = (e_vec/e_mag)*cos(phi_s) + (cross(h_hat, e_vec)/norm(cross(h_hat,e_vec)))*sin(phi_s)
e_hat = e_vec / e_mag;
h_cross_e_unit = cross(h_hat, e_vec) / norm(cross(h_hat, e_vec));
S_hat = e_hat * cos(phi_s) + h_cross_e_unit * sin(phi_s);

% T_hat is the "transverse" or "target" axis in the B-plane
% It's often aligned with the ecliptic or equatorial plane.
S_cross_K = cross(S_hat, K_hat);
T_hat = S_cross_K / norm(S_cross_K);

% R_hat completes the right-handed coordinate system (S, T, R)
R_hat = cross(S_hat, T_hat);

% --- Step 4: Calculate the B-Vector and its Components ---

% The B-vector lies in the B-plane, is perpendicular to S_hat, and has magnitude b
B_vec = b * cross(S_hat, h_hat);

% Project the B-vector onto the T and R axes to get the final components
B_T = dot(B_vec, T_hat);
B_R = dot(B_vec, R_hat);

end