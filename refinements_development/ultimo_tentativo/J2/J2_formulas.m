% MATLAB script to compute and use the J2 Gravity Gradient Matrix for Saturn
clear all
close all
clc
%% 1. Define Symbolic Variables
syms x y z mu J2 R_ref real

%% 2. Define J2 Acceleration and Position Vectors
r_vec = [x; y; z];
r_norm = sqrt(x^2 + y^2 + z^2);

% Common factor from the J2 acceleration formula
common_factor = - (3 * mu * J2 * R_ref^2) / (2 * r_norm^5);

% J2 acceleration vector components
z_squared_over_r_squared = z^2 / r_norm^2;
ax = common_factor * x * (1 - 5 * z_squared_over_r_squared);
ay = common_factor * y * (1 - 5 * z_squared_over_r_squared);
az = common_factor * z * (3 - 5 * z_squared_over_r_squared);

a_J2 = [ax; ay; az];

%% 3. Compute the Symbolic Gravity Gradient Matrix (Jacobian)
fprintf('Computing the symbolic Jacobian matrix...\n');
G_J2_symbolic = jacobian(a_J2, r_vec);
fprintf('Computation complete.\n');

% Display a portion of the symbolic matrix to see the result
% Note: The full expression is very long
fprintf('Example partial derivative (d(ax)/dx):\n');
disp(G_J2_symbolic(1,1));

%% 4. Convert to a Numerical MATLAB Function
fprintf('Converting symbolic matrix to a numerical function handle...\n');
G_J2_function = matlabFunction(G_J2_symbolic, 'vars', {[x; y; z], mu, J2, R_ref});
fprintf('Function handle created.\n\n');

%% 5. Example Usage with Saturn's Constants
% Define Saturn's parameters
mu_Saturn = 3.7931187e7; % km^3/s^2
J2_Saturn = 1.629071e-2;  % Nondimensional
R_ref_Saturn = 60330;      % km (Equatorial radius)

% Example spacecraft position vector (in km)
spacecraft_position = [70000; 25000; -15000];

% Calculate the numerical gravity gradient matrix at this position
G_J2_numerical = G_J2_function(spacecraft_position, mu_Saturn, J2_Saturn, R_ref_Saturn);

fprintf('Numerical J2 Gravity Gradient Matrix at position [%g, %g, %g] km:\n', ...
        spacecraft_position(1), spacecraft_position(2), spacecraft_position(3));
disp(G_J2_numerical);