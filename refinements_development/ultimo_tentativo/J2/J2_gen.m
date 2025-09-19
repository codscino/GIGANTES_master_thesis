% --- generate_J2_gradient_function.m ---
%
% This script uses the Symbolic Math Toolbox to derive the analytical J2
% gravity gradient and saves it as an optimized .m file for fast execution.
% Run this script only ONCE to create the function file.

clear; clc;
fprintf('Starting symbolic derivation of the J2 gravity gradient...\n');

%% 1. Define Symbolic Variables
% Define inputs to the future function: position vector and parameters
syms x y z real
syms mu J2 R_ref real

%% 2. Formulate the J2 Acceleration Vector
r_vec = [x; y; z];
r_norm = sqrt(x^2 + y^2 + z^2);

common_factor = - (3 * mu * J2 * R_ref^2) / (2 * r_norm^5);
z_squared_over_r_squared = z^2 / r_norm^2;

ax = common_factor * x * (1 - 5 * z_squared_over_r_squared);
ay = common_factor * y * (1 - 5 * z_squared_over_r_squared);
az = common_factor * z * (3 - 5 * z_squared_over_r_squared);

a_J2 = [ax; ay; az];

%% 3. Compute the Symbolic Gravity Gradient Matrix (Jacobian)
G_J2_symbolic = jacobian(a_J2, r_vec);

% For numerical stability and efficiency, it's often better to simplify
% the expressions. This can take a moment.
fprintf('Simplifying symbolic expressions...\n');
G_J2_symbolic = simplify(G_J2_symbolic);

%% 4. Generate the Optimized MATLAB Function File
% This is the key step. We will create a new file named
% 'computeJ2GravityGradient_generated.m'.
output_filename = 'computeJ2GravityGradient_generated.m';

fprintf('Generating the MATLAB function file: %s\n', output_filename);
matlabFunction(G_J2_symbolic, ...
               'File', output_filename, ...
               'Vars', {[x; y; z], mu, J2, R_ref}, ...
               'Outputs', {'G_J2'});

fprintf('Script finished. The function %s has been created.\n', output_filename);