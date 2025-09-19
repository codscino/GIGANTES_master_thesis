function [tt, yy, te, ye, ie] = propagateKeplerODE2(r0, v0, timevector, mu, event_func_handle)
    % propagateKeplerODE2: Propagates a two-body orbit using MATLAB's ODE45.
    %
    % This function integrates the Keplerian equations of motion for a given
    % initial state (position and velocity) over a specified time vector. It
    % can optionally detect specific events during the integration, such as
    % crossing a certain altitude or sphere of influence.
    %
    % SYNTAX:
    %   [tt, yy] = propagateKeplerODE2(r0, v0, timevector, mu)
    %   [tt, yy, te, ye, ie] = propagateKeplerODE2(r0, v0, timevector, mu, event_func_handle)
    %
    % INPUTS:
    %   r0          - [3x1] Initial position vector of the spacecraft (km).
    %   v0          - [3x1] Initial velocity vector of the spacecraft (km/s).
    %   timevector  - [Nx1] Vector of time points at which to report the
    %                 solution (s). The integration starts at timevector(1)
    %                 and ends at timevector(end).
    %   mu          - [1x1] Gravitational parameter of the central body
    %                 (km^3/s^2).
    %   event_func_handle - (Optional) [function handle] A handle to an ODE
    %                 events function, created using the same structure as
    %                 MATLAB's odeset 'Events' property. If not provided or
    %                 left empty, no event detection will be performed.
    %
    % OUTPUTS:
    %   tt          - [Mx1] Output time vector corresponding to the rows of yy.
    %                 If the integration is terminated early by an event, M
    %                 will be less than N.
    %   yy          - [Mx6] Solution matrix, where each row corresponds to a
    %                 time in tt. The columns are [rx, ry, rz, vx, vy, vz].
    %   te          - [Px1] Vector of times at which events occurred. Empty if
    %                 no events were detected.
    %   ye          - [Px6] Matrix of spacecraft states [r; v] at the event
    %                 times in te.
    %   ie          - [Px1] Vector of indices indicating which event was
    %                 triggered at the corresponding time in te.
    %

    % By default, nargin is 5 if the argument is passed, even if it's empty.
    % So we check if it was not passed or if it is empty.
    if nargin < 5 || isempty(event_func_handle)
        has_events = false;
    else
        has_events = true;
    end

    F = @(t, x) keplerEOM(t, x, mu);

    % --- Set base options ---
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

    % --- Call the solver conditionally ---
    if has_events
        % If events are provided, add them to options and call ode45 requesting all outputs.
        options = odeset(options, 'Events', event_func_handle);
        [tt, yy, te, ye, ie] = ode45(F, timevector, [r0; v0], options);
    else
        % If no events, call ode45 requesting only t and y.
        [tt, yy] = ode45(F, timevector, [r0; v0], options);
        % Assign empty arrays to the event outputs for consistent function signature.
        te = [];
        ye = [];
        ie = [];
    end
end

% --- Nested function for Keplerian Equations of Motion ---
function dxdt = keplerEOM(~, x, mu)
    r = x(1:3);
    a = -mu * r / norm(r)^3;
    dxdt = [x(4:6); a];
end