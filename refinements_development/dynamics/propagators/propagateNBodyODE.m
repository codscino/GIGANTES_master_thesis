function [tt, yy] = propagateNBodyODE(rvec, vvec, timevector, muCB, mu_TBs, ephem_func, useJ2)
% Propagate a point-mass in an n-body system, with an optional J2 perturbation.
%
%   Inputs
%   ======
%   ... (rvec, vvec, timevector, muCB are the same) ...
%   useJ2       (logical) Flag to enable/disable J2 perturbation (true/false).
%   J2_CB       (double)  J2 coefficient of the central body.
%   R_eq_CB     (double)  Equatorial radius of the central body (km).
%   mu_TBs      (vector)  Gravitational parameters for M target bodies.
%   ephem_func  (handle)  Function handle to get target body positions.

if nargin < 7
    useJ2 = false;
end

% --- Pass all parameters to the Equations of Motion (EOM) ---
F = @(t, x) nBodyEOM(t, x, muCB, mu_TBs, ephem_func, useJ2);

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[tt, yy] = ode45(F, timevector, [rvec; vvec], options);

end

% --- Nested function for the Equations of Motion ---
function dxdt = nBodyEOM(t, x, muCB, mu_TBs, ephem_func, useJ2)
    % Unpack the state vector
    r_sc = x(1:3); % Spacecraft position vector relative to Central Body

    % 1. Acceleration from the Central Body (as a point mass)
    a_CB = -muCB * r_sc / (norm(r_sc)^3);

    % 2. Acceleration from J2 Perturbation (if enabled)
    a_J2 = [0; 0; 0]; % Initialize to zero
    if useJ2
        a_J2 = calculate_J2_acceleration(r_sc);
    end

    % 3. Get positions of all perturbing Target Bodies (TBs) at time t
    r_TBs = ephem_func(t);
    
    % Calculate the third-body perturbation sum
    a_pert_sum = [0; 0; 0]; 
    for i = 1:length(mu_TBs)
        mu_tb = mu_TBs(i);
        r_tb  = r_TBs(:, i);
        rtb2sc = r_sc - r_tb;
        
        % The direct term, a_direct = -mu_tb * rtb2sc / |rtb2sc|^3
        a_direct_i = - mu_tb * (rtb2sc / norm(rtb2sc)^3);
        % The indirect term, a_indirect = +mu_tb * r_tb / |r_tb|^3
        a_indirect_i =  + mu_tb*(r_tb / norm(r_tb)^3);
        % The total perturbation is (a_direct - a_indirect).
        pert_term_i =  a_direct_i - a_indirect_i;
        
        a_pert_sum = a_pert_sum + pert_term_i;
    end
    
    % 4. Total Acceleration is the sum of all forces
    a_total = a_CB + a_J2 + a_pert_sum;
    
    % 5. Assemble the state derivative vector
    dxdt = [x(4:6); a_total];
end