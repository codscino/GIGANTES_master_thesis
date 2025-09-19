% This script simulates the motion of a bouncing ball using
% the ode45 solver and an event function to detect bounces.

% --- 1. Set Up Figure ---
clear; clc; close all; 

figure;      % Create a new figure window for the plot
hold on;     % Hold the plot to draw multiple bounce trajectories
grid on;
title('Bouncing Ball Simulation');
xlabel('Time (s)');
ylabel('Height (m)');

% --- 2. Define Initial Conditions and Simulation Parameters ---
tstart = 0;          % Start time of the first bounce
tfinal = 30;         % Maximum simulation time
y0 = [0; 20];        % Initial conditions: [position; velocity] -> [0m; 20m/s]
max_bounces = 10;    % Set a limit for the number of bounces

% --- 3. Set ODE Solver Options ---
options = odeset('Events', @bounceEvents);

% --- 4. Main Simulation Loop ---
% This loop runs for each bounce, restarting the solver with new initial conditions.
for i = 1:max_bounces
    % Call the ode45 solver. It will integrate from tstart to tfinal,
    % or until the event defined in 'bounceEvents' occurs.
    % te, ye, ie will contain information about the event if one is found.
    time_step = 0.01;
    tspan = tstart:time_step:tfinal;
    [t, y, te, ye, ie] = ode45(@ballODE, tspan, y0, options);

    % Plot the trajectory of the current bounce
    plot(t, y(:,1), 'b-'); % Plot position (first column of y) vs. time
    drawnow; % Update the plot window immediately

    % Check if an event (a bounce) was detected. If not, the ball
    % didn't return to the ground, so we can stop the simulation.
    if isempty(te)
        disp('Ball did not return to the ground. Simulation ended.');
        break;
    end

    % --- 5. Update Conditions for the Next Bounce ---
    % Set the start time for the next integration period to the time of the last bounce.
    tstart = te;

    % Set the new initial conditions for the next bounce.
    % The initial position is the ground (0).
    % The initial velocity is the velocity at impact (ye(2)), reversed
    % and damped by a factor of 0.9 to simulate energy loss.
    y0(1) = 0;
    y0(2) = -0.6 * ye(2);

    % Optional: Add a small check to stop if the bounces are too small
    if abs(y0(2)) < 1e-2
        disp('Ball has come to rest. Simulation ended.');
        break;
    end
end

hold off; % Release the plot hold
disp('Simulation complete.');


%% Functions

function dydt = ballODE(t, y)
    % This function defines the differential equations for the ball's motion.
    % y(1) is position, y(2) is velocity.
    % The derivative of position is velocity.
    % The derivative of velocity is acceleration (gravity = -9.8).
    dydt = [y(2); -9.8];
end

function [value, isterminal, direction] = bounceEvents(t, y)
    % This function detects the event: the ball hitting the ground.
    % An event occurs when 'value' is zero.
    value = y(1);      % We want to detect when the height, y(1), is zero.

    % 'isterminal' = 1 tells the solver to stop the integration when the event occurs.
    isterminal = 1;

    % 'direction' = -1 tells the solver to only detect the event if the 'value'
    % function is decreasing (i.e., the ball is falling).
    direction = -1;
end