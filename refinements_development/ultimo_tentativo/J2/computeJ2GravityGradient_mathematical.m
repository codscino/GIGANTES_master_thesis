function G_J2 = computeJ2GravityGradient_mathematical(r, mu, J2, R_ref)
%computeJ2GravityGradient_mathematical - J2 gravity gradient using mathematical formulation
%
% INPUTS:
%   r     - [3x1] Position vector (km)
%   mu    - Gravitational parameter (km^3/s^2)
%   J2    - J2 coefficient (dimensionless)
%   R_ref - Reference radius (km)
%
% OUTPUT:
%   G_J2  - [3x3] J2 gravity gradient matrix
%
% This implementation uses the mathematical formulation with sign convention
% consistent with G = ∂a/∂r for direct addition to other gradient terms

% Extract position components
x = r(1);
y = r(2);
z = r(3);

% Compute radius and powers
r_mag = norm(r);
r2 = r_mag^2;
r5 = r_mag^5;

% Common prefactor (negative to match the required sign convention)
prefactor = -(3/2) * mu * J2 * R_ref^2 / r5;

% Compute normalized position ratios
x2_r2 = x^2 / r2;
y2_r2 = y^2 / r2;
z2_r2 = z^2 / r2;

% Matrix components using the formulation that matches your code's convention
% These expressions are derived from -G' where G' is the standard mathematical form

% Diagonal terms
G_xx = -1 + 7*x2_r2 + 5*z2_r2 - 25*x2_r2*z2_r2;
G_yy = -1 + 7*y2_r2 + 5*z2_r2 - 25*y2_r2*z2_r2;
G_zz = -3 + 15*z2_r2 - 25*z2_r2*z2_r2;

% Off-diagonal terms (symmetric)
G_xy = 5*x*y/r2 * (1 - 7*z2_r2);
G_xz = 5*x*z/r2 * (3 - 7*z2_r2);
G_yz = 5*y*z/r2 * (3 - 7*z2_r2);

% Assemble the symmetric matrix
G_J2 = prefactor * [
    G_xx,  G_xy,  G_xz;
    G_xy,  G_yy,  G_yz;
    G_xz,  G_yz,  G_zz
];

end