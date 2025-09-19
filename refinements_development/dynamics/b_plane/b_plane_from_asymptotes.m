function [B_vec, r_p, B_R, B_T] = b_plane_from_asymptotes(v_inf_in, v_inf_out, mu)
% Calculates B-plane parameters from asymptote vectors.
%
% From algorithm 79 of Vallado
%
% Inputs:
%   v_inf_in    (3x1) Incoming hyperbolic excess velocity vector [km/s]
%   v_inf_out   (3x1) Outgoing hyperbolic excess velocity vector [km/s]
%   mu          (1x1) Gravitational parameter of the central body [km^3/s^2]
%
% Outputs:
%   B_vec       (3x1) The B-vector [km]
%   r_p         (1x1) Radius of periapsis (closest approach) [km]
%   B_R         (1x1) Component of the B-vector along the R-axis [km] (optional)
%   B_T         (1x1) Component of the B-vector along the T-axis [km] (optional)
%
% Note:
% For an unpowered gravity assist, the magnitudes of v_inf_in and
% v_inf_out should be equal. The function does not enforce this, allowing
% for analysis of powered flybys as well. The incoming velocity magnitude
% is used as the reference for calculations.
%

% --- Input Validation ---
% Ensure inputs are column vectors for consistent matrix operations
v_inf_in = v_inf_in(:);
v_inf_out = v_inf_out(:);

% --- Step 1: Define the Orbit and B-Plane Reference Frames ---

% S_hat is the unit vector along the incoming asymptote
v_inf_in_mag = norm(v_inf_in);
S_hat = v_inf_in / v_inf_in_mag;

% h_hat is the unit vector normal to the orbit plane
h_vec = cross(v_inf_in, v_inf_out);
h_hat = h_vec / norm(h_vec);

% B_hat is the unit vector in the direction of the B-vector
B_hat = cross(S_hat, h_hat);

% K_hat is the fundamental reference direction (e.g., Z-axis of ECI frame)
K_hat = [0; 0; 1];

% T_hat is the "transverse" or "target" axis in the B-plane.
% Note: The algorithm lists T = S x K, but T_hat must be a unit vector.
T_vec = cross(S_hat, K_hat);
T_hat = T_vec / norm(T_vec);

% R_hat completes the right-handed coordinate system (S, T, R)
R_hat = cross(S_hat, T_hat);


% --- Step 2: Calculate Flyby Geometry (Turn Angle and Periapsis) ---

% Calculate the turn angle (phi) between the asymptotes
cos_phi = dot(v_inf_in, v_inf_out) / (v_inf_in_mag * norm(v_inf_out));
% Clamp the value to avoid numerical errors with acos for angles near 0 or 180
cos_phi = max(-1, min(1, cos_phi)); 
phi = acos(cos_phi);

% Calculate the radius of periapsis (r_p)
v_inf_in_mag_sq = v_inf_in_mag^2;
r_p_term_denom = cos((pi - phi)/2);
r_p = (mu / v_inf_in_mag_sq) * ( (1 / r_p_term_denom) - 1 );


% --- Step 3: Calculate the B-Vector ---

% First, calculate the magnitude of the B-vector (B)
% This formula comes from combining B=a*sqrt(e^2-1) and e = 1 + (rp*v_inf^2)/mu
term_inside_sqrt = (1 + (v_inf_in_mag_sq * r_p) / mu)^2 - 1;
B_mag = (mu / v_inf_in_mag_sq) * sqrt(term_inside_sqrt);

% Construct the full B-vector
B_vec = B_mag * B_hat;

% --- Step 4: Calculate Optional B-Plane Components ---
% Project the B-vector onto the T and R axes
if nargout > 2
    B_T = dot(B_vec, T_hat);
    B_R = dot(B_vec, R_hat);
end

end