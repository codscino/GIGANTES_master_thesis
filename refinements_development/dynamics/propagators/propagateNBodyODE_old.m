function [tt, yy] = propagateNBodyODE(rvec, vvec, timevector, muCB, mu_TBs, ephem_func)
% Propagate a point-mass in an n-body system.
%
%   Inputs
%   ======
%   rvec        3x1 initial position of spacecraft w.r.t. central body (km)
%   vvec        3x1 initial velocity of spacecraft w.r.t. central body (km/s)
%   timevector  Nx1 vector of time points (s)
%   muCB        Scalar gravitational parameter of central body (km^3/s^2)
%   mu_TBs      Mx1 vector of gravitational parameters for M target bodies
%   ephem_func  Function handle that takes time 't' and returns a 3xM
%               matrix of target body positions w.r.t. the central body.
%               Example: r_moons = ephem_func(t);
%
%   Outputs
%   =======
%   tt          Nx1 output times (s)
%   yy          Nx6 matrix of state history [x, y, z, vx, vy, vz]

    % anonymous function for the equations of motion
    F = @(t, x) nBodyEOM(t, x, muCB, mu_TBs, ephem_func);
    
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); % Higher precision for n-body
    [tt, yy] = ode45(F, timevector, [rvec; vvec], options);
    
    end
    
    % --- Nested function for the Equations of Motion ---
    function dxdt = nBodyEOM(t, x, muCB, mu_TBs, ephem_func)
        % Unpack the state vector
        r_sc = x(1:3); % Spacecraft position vector relative to Central Body
        
        % 1. Acceleration from the Central Body
        a_CB = -muCB * r_sc / (norm(r_sc)^3);
        
        % 2. Get positions of all perturbing Target Bodies (TBs) at time t
        r_TBs = ephem_func(t); % This returns a 3xM matrix
        
        % 3. Calculate the perturbation term to be subtracted
        % This loop calculates a_pert = SUM { mu_tb * (rtb2sc/|...|^3 + r_tb/|...|^3) }
        a_pert_sum = [0; 0; 0]; 
        for i = 1:length(mu_TBs)
            mu_tb = mu_TBs(i);
            r_tb  = r_TBs(:, i); % Position of i-th target body relative to CB
        
            % Vector from Target Body to Spacecraft
            rtb2sc = r_sc - r_tb;
        
            % The direct term, a_direct = -mu_tb * rtb2sc / |rtb2sc|^3
            a_direct_i = - mu_tb * (rtb2sc / norm(rtb2sc)^3);
            % The indirect term, a_indirect = +mu_tb * r_tb / |r_tb|^3
            a_indirect_i =  + mu_tb*(r_tb / norm(r_tb)^3);
            % The total perturbation is (a_direct - a_indirect).
            pert_term_i =  a_direct_i - a_indirect_i;
            
            a_pert_sum = a_pert_sum + pert_term_i;
        end
        
        % 4. Total Acceleration a_total = a_CB + a_pert_sum
        a_total = a_CB + a_pert_sum;
        
        % 5. Assemble the state derivative vector
        dxdt = [x(4:6); a_total];

end