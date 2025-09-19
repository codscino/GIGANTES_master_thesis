function a_out = calculate_zonal_perturbations(r_vec, type)
% Computes accelerations from J2 to J6 for Saturn
%
%
% INPUTS:
%   r_vec   (3x1 double) The position vector of the spacecraft relative to
%           the central body [x; y; z] in km.
%   type    (1xN int array) An array of integers specifying which zonal
%           harmonics to compute. e.g., [2, 3] for J2 and J3.
%           Valid entries are 2, 3, 4, 5, 6.
%
% OUTPUTS:
%   a_out   (3xN double) A matrix where each column corresponds to the
%           perturbing acceleration vector for the requested J-term.
%           For example, if type = [2, 4], a_out will be a 3x2 matrix
%           with [a_J2, a_J4].
%
%           To get the total combined acceleration, use:
%           a_total = sum(a_out, 2);
%

% --- Constants for Saturn ---
mu = getAstroConstants('Saturn', 'Mu');
J2 = 16290.573/10^6; % 
r_eq = 60268; % equatorial Saturn radius(values from Spice)

% Zonal Harmonic Coefficients (unnormalized) source for J2-J6:
% https://atmos.nmsu.edu/data_and_services/atmospheres_data/Cassini/logs/Table%201Measured%20gravity%20harmonic%20coefficients.pdf
J_coeffs = [
    0;          % J0 (not used)
    0;          % J1 (not used)
    16290.573  % J2 16290.573
    0.059;      % J3
    -935.314;   % J4
    -0.224;     % J5
    86.340;     % J6
]/10^6; % normalize them

% --- Initial Calculations ---
x = r_vec(1);
y = r_vec(2);
z = r_vec(3);

r_mag = norm(r_vec);
r_mag2 = r_mag^2;
r_eq_over_r = r_eq / r_mag;

% Pre-calculate powers of z/r for efficiency
s = z / r_mag; % This is sin(latitude)
s2 = s^2;
s3 = s * s2;
s4 = s * s3;
s5 = s * s4;
s6 = s * s5;

% --- Loop and Calculate Accelerations ---

% Initialize output matrix
num_requests = length(type);
a_out = zeros(3, num_requests);

% Common scalar prefix for all terms
common_prefix = -mu / r_mag2;

for k = 1:num_requests
    n = abs(type(k)); % Use absolute value to be robust (e.g., handles -2)
    
    a_J = zeros(3,1); % Initialize current acceleration vector
    
    % Select the correct formula based on the J-term number 'n'
    switch n
        case 2 % J2 
            C2 = 1.5 * J_coeffs(n+1) * r_eq_over_r^2;
            
            gamma = 1 - 5*s2; % Common factor for x and y components
            ax = C2 * gamma * (x/r_mag);
            ay = C2 * gamma * (y/r_mag);
            az = C2 * (3 - 5*s2) * s;
            
            a_J = [ax; ay; az];

        case 3 % J3
            C3 = 0.5 * J_coeffs(n+1) * r_eq_over_r^3;
            ax = C3 * 5 * (3 - 7*s2) * s * (x/r_mag);
            ay = C3 * 5 * (3 - 7*s2) * s * (y/r_mag);
            az = C3 * (35*s4 - 30*s2 + 3) * s;  % Added the s factor
            a_J = [ax; ay; az];

        case 4 % J4 from Eq (10.66)
            C4 = 1.25 * J_coeffs(n+1) * r_eq_over_r^4; % 5/4, not 5/8 as in image
                                                      % 5/8 is a common typo, Vallado uses 5/4
            
            gamma = 3 - 42*s2 + 63*s4;
            ax = C4 * gamma * (x/r_mag);
            ay = C4 * gamma * (y/r_mag);
            az = C4 * (-15 + 70*s2 - 63*s4) * s;

            a_J = [ax; ay; az];
            
        case 5 % J5 from Eq (10.67)
            C5 = (1/8) * J_coeffs(n+1) * r_eq_over_r^5; % Use Vallado's constant
            ax = C5 * 3*s * (35 - 210*s2 + 231*s4) * (x/r_mag);
            ay = C5 * 3*s * (35 - 210*s2 + 231*s4) * (y/r_mag);
            az = C5 * (693*s6 - 945*s4 + 315*s2 - 15) * s;
            a_J = [ax; ay; az];
     

        case 6 % J6 from Eq (10.68)
            C6 = (1/16) * J_coeffs(n+1) * r_eq_over_r^6;
            
            gamma = 35 - 945*s2 + 3465*s4 - 3003*s6;
            ax = C6 * gamma * (x/r_mag);
            ay = C6 * gamma * (y/r_mag);
            az = C6 * (3003*s6 - 4851*s4 + 2205*s2 - 315) * s;
            
            a_J = [ax; ay; az];
            
        otherwise
            warning('Zonal harmonic J%d is not implemented. Skipping.', n);
            % The column for this request will remain zeros
    end
    
    % The full acceleration is the calculated vector multiplied by the common prefix
    a_out(:, k) = common_prefix * a_J;
    
end

end