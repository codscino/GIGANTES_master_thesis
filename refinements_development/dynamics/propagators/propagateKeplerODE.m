function [tt, yy] = propagateKeplerODE(rvec, vvec, timevector, muCB)
%PROPAGATEKEPLERODE Propagate a point‐mass under Keplerian gravity.
%
%   Inputs
%   ======
%   rvec        3×1 column vector of initial position [x; y; z] (km)
%   vvec        3×1 column vector of initial velocity [vx; vy; vz] (km/s)
%   timevector  N×1 column vector of time points at which to compute the state (s)
%   muCB         scalar gravitational parameter of central body (km^3/s^2)
%
%   Outputs
%   =======
%   tt          N×1 column vector of output times (s), same as timevector
%   yy          N×6 matrix of state history, where each row is
%                   [x, y, z, vx, vy, vz]
%
%   Example
%   -------
%       % Propagate from t=0 to 3600 s in 60-s steps:
%       t = (0:60:3600).';
%       [t_out, y_out] = propagateKeplerODE([7000;0;0],[0;7.5;0], t, 398600);
%       % y_out(:,1:3) are positions, y_out(:,4:6) are velocities
%
%Equation of motion ¨r + mu*r/|r|^3 = 0
F = @(t,x) [ ...
    x(4);                    % Vx
    x(5);                    % Vy
    x(6);                    % Vz
   -muCB*x(1)/(norm(x(1:3))^3);  % ax
   -muCB*x(2)/(norm(x(1:3))^3);  % ay
   -muCB*x(3)/(norm(x(1:3))^3)]; % az

options = odeset('RelTol',1e-6,'AbsTol',1e-7,'Refine',50);
[tt,yy] = ode45(F, timevector, [rvec; vvec], options);

end