function [ra_linked_conics, ra_kepler, ra_nbody, f_in] = calculateApoapsisComparison(pars, spiceParam, t0_mjd, perturbingBodyNaifIDs,  propDirection, useJ2)
% Propagates an orbit and compares apoapsis results.
%
% INPUTS:
%   pars                (struct)  A structure containing simulation parameters.
%                                 Must contain the following fields:
%                                 .INPUTS.idCentral (int): Planet-centric ID for the central body (e.g., 6 for Saturn).
%                                 .INPUTS.idMoon (int): Moon-centric ID for the flyby body (e.g., 1 for Enceladus).
%                                 .INPUTS.V_inf (double): Hyperbolic excess velocity for the flyby [km/s].
%                                 .INPUTS.alpha (double): Crank angle [rad].
%                                 .INPUTS.k (double): Aiming parameter 'k'.
%
%   kernels             (cell)    A cell array of strings with the names of the
%                                 SPICE kernels required for ephemeris data.
%                                 e.g., {'sat441.bsp', 'naif0012.tls', 'pck00010.tpc', 'de440.bsp'}
%
%   spiceParam          (struct)  A structure with SPICE ephemeris parameters.
%                                 .frame (char): The reference frame (e.g., 'J2000').
%                                 .abcorr (char): Aberration correction (e.g., 'NONE').
%                                 .observer (char): The NAIF ID of the observer (e.g., '699' for Saturn).
%
%   t0_mjd              (double)  The initial epoch for the simulation, given
%                                 as Modified Julian Date 2000.
%
%   perturbingBodyNaifIDs (vector) A row vector of NAIF IDs for the third-body
%                                  perturbers in the N-body simulation.
%                                  e.g., [601, 602, 603, 604, 605, 606, 10] for
%                                  major Saturn moons and the Sun.
%
%   propDirection         (string) Either 'forward' or 'backward'     
%
% OUTPUTS:
%   a_linked_conics     (double)  The calculated apoapsis radius from the
%                                 analytical linked-conics model [km].
%   a_kepler            (double)  The calculated apoapsis radius from the
%                                 Keplerian (2-body) propagation [km].
%   a_nbody             (double)  The calculated apoapsis radius from the
%                                 N-body propagation [km].

if nargin < 6
    useJ2 = false;
end

if nargin < 5
    propDirection = 'forward'; % Default to forward propagation
end

%% ========================================================================
%  1. DEFINE COMMON PARAMETERS AND INITIAL STATE
%  ========================================================================

% --- Load Constants based on Inputs ---
[pars.Planet.mu, ~, ~, ~] = planetConstants(pars.INPUTS.idCentral);
[~, pars.Moon.mu, ~, ~] = satMoonsConstants(pars.INPUTS.idMoon);

% --- Generate Initial State from Flyby Parameters ---
pars.INPUTS.epoch0 = 0; % Relative epoch for flyby calculation

[~, r0_sc, v0_sc, ~] = vinfAlphaCrank_to_VinfCART(pars.INPUTS.V_inf, pars.INPUTS.alpha,...
    pars.INPUTS.k, pars.INPUTS.epoch0, pars.INPUTS.idMoon, pars.INPUTS.idCentral);

%% ========================================================================
%  2. ANALYSIS A: IDEALIZED LINKED CONICS (ANALYTICAL)
%  ========================================================================

mu_central_body = pars.Planet.mu;

initial_state = [r0_sc, v0_sc];
[kep1] = car2kep(initial_state, mu_central_body);
f_in = rad2deg(kep1(6));

% ra = a*(1+e)
a = kep1(1);
e = kep1(2);
ra_linked_conics = a * (1 + e);

% Define propagation duration based on the analytical orbit period
lc_period_sec = 2*pi * sqrt(kep1(1)^3 / mu_central_body); 
duration_sec = lc_period_sec * 1.2; % Propagate for 1.2 orbits

% Create time vector based on propagation direction
if strcmpi(propDirection, 'forward')
    time_vector = linspace(0, duration_sec, 1000)';
elseif strcmpi(propDirection, 'backward')
    time_vector = linspace(0, -duration_sec, 1000)'; % just a minus on duration_sec
else
    error('Invalid propDirection specified. Use "forward" or "backward".');
end

%% ========================================================================
%  3. ANALYSIS B: KEPLERIAN (2-BODY) PROPAGATION
%  ========================================================================

[~, kep_state] = propagateKeplerODE(r0_sc', v0_sc', time_vector, mu_central_body);

kep_r_norms = vecnorm(kep_state(:, 1:3), 2, 2);
[ra_kepler, ~] = max(kep_r_norms);


%% ========================================================================
%  4. ANALYSIS C: N-BODY PROPAGATION
%  ========================================================================

% --- Dynamically build the list of gravitational parameters (mu) ---
mu_TBs = zeros(1, length(perturbingBodyNaifIDs));
for i = 1:length(perturbingBodyNaifIDs)
    id = perturbingBodyNaifIDs(i);
    switch id
        case 10 % Sun
            mu_TBs(i) = getAstroConstants('Sun', 'Mu');
        case 599 % Jupiter
            mu_TBs(i) = getAstroConstants('Jupiter', 'Mu');
        case 601 % Mimas
            [~, mu, ~, ~] = satMoonsConstants(0); mu_TBs(i) = mu;
        case 602 % Enceladus
            [~, mu, ~, ~] = satMoonsConstants(1); mu_TBs(i) = mu;
        case 603 % Tethys
            [~, mu, ~, ~] = satMoonsConstants(2); mu_TBs(i) = mu;
        case 604 % Dione
            [~, mu, ~, ~] = satMoonsConstants(3); mu_TBs(i) = mu;
        case 605 % Rhea
            [~, mu, ~, ~] = satMoonsConstants(4); mu_TBs(i) = mu;
        case 606 % Titan
            [~, mu, ~, ~] = satMoonsConstants(5); mu_TBs(i) = mu;
        case 607 % Hyperion
            mu_TBs(i) = 0.374; % From script
        case 608 % Iapetus
            mu_TBs(i) = 120.4; % From script
        otherwise
            error('calculateApoapsisComparison:UnknownNaifID', ...
                  'The NAIF ID %d is not recognized for mu calculation.', id);
    end
end
      
%--- Setup Ephemeris Handle and Propagate ---
ephem_handle = @(t_sec) get_all_moon_positions(t_sec, t0_mjd, perturbingBodyNaifIDs, spiceParam);

[~, nbody_state] = propagateNBodyODE(r0_sc', v0_sc', time_vector, ...
                                   mu_central_body, mu_TBs, ephem_handle, useJ2);

nbody_r_norms = vecnorm(nbody_state(:, 1:3), 2, 2);
[ra_nbody, ~] = max(nbody_r_norms);



end