function [a_lc, a_nbody, f_in] = lc_vs_nbody_sma2(pars, spiceParam, t0_mjd, perturbingBodyNaifIDs,  propDirection)
% Propagates an orbit and compares semi major axis results.
%
% INPUTS:
%   pars                (struct)  A structure containing simulation parameters.
%                                 Must contain the following fields:
%                                 .INPUTS.idCentral (int): Planet-centric ID for the central body (e.g., 6 for Saturn).
%                                 .INPUTS.idMoon (int): Moon-centric ID for the flyby body (e.g., 1 for Enceladus).
%                                 .INPUTS.V_inf (double): Hyperbolic excess velocity for the flyby [km/s].
%                                 .INPUTS.alpha (double): Pump angle [rad].
%                                 .INPUTS.k (double): Crank angle [rad].
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
%   perturbingBodyNaifIDs (vector) A row vector of NAIF IDs for third-body
%                                  perturbers and special perturbation flags.
%                                  - Body NAIF IDs (e.g., 10 for Sun, 606 for Titan).
%                                  - Special IDs:
%                                    -2: Include J2 perturbation.
%                                    -3: Include J3 perturbation.
%                                    -4: Include J4 perturbation.
%                                  e.g., [606, 10, -2, -4] for Titan, Sun, J2, and J4.
%
%   propDirection         (string) Either 'forward' or 'backward'     
%
% OUTPUTS:
%   a_lc                (double)  The calculated sma radius from the
%                                 analytical linked-conics model [km].
%
%   a_nbody             (double)  The calculated sma from the
%                                 N-body propagation [km].
%


if nargin < 5
    propDirection = 'forward';
end

%% ========================================================================
%  1. PARSE INPUTS AND SETUP
%  ========================================================================

% --- Filter the input IDs into actual bodies and special perturbations ---
actualBodyNaifIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs >= 0);
specialPerturbationIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs < 0);

% --- Load Constants & Generate Initial State ---
[pars.Planet.mu, ~, ~, ~] = planetConstants(pars.INPUTS.idCentral);
mu_central_body = pars.Planet.mu;

[~, r0_sc, v0_sc, ~] = vinfAlphaCrank_to_VinfCARTClaudio(pars.INPUTS.V_inf, pars.INPUTS.alpha,...
    pars.INPUTS.k, t0_mjd, pars);

%% ========================================================================
%  2. ANALYSIS A: IDEALIZED LINKED CONICS (ANALYTICAL)
%  ========================================================================
initial_state = [r0_sc, v0_sc];
[kep1] = car2kep(initial_state, mu_central_body);
f_in = rad2deg(kep1(6));
a_lc = kep1(1);

lc_period_sec = 2*pi * sqrt(a_lc^3 / mu_central_body); 
duration_sec = lc_period_sec * 1.2;

if strcmpi(propDirection, 'forward')
    time_vector = linspace(0, duration_sec, 1000)';
elseif strcmpi(propDirection, 'backward')
    time_vector = linspace(0, -duration_sec, 1000)';
else
    error('Invalid propDirection specified. Use "forward" or "backward".');
end

%% ========================================================================
%  3. ANALYSIS C: N-BODY PROPAGATION
%  ========================================================================

% --- Dynamically build the list of gravitational parameters (mu) for actual bodies ---
mu_TBs = zeros(1, length(actualBodyNaifIDs));
for i = 1:length(actualBodyNaifIDs)
    id = actualBodyNaifIDs(i);
    switch id
        case 10,  mu_TBs(i) = getAstroConstants('Sun', 'Mu');
        case 5, mu_TBs(i) = getAstroConstants('Jupiter', 'Mu');
        case 601, [~, mu, ~, ~] = satMoonsConstants(0); mu_TBs(i) = mu; % Mimas
        case 602, [~, mu, ~, ~] = satMoonsConstants(1); mu_TBs(i) = mu; % Enceladus
        case 603, [~, mu, ~, ~] = satMoonsConstants(2); mu_TBs(i) = mu; % Tethys
        case 604, [~, mu, ~, ~] = satMoonsConstants(3); mu_TBs(i) = mu; % Dione
        case 605, [~, mu, ~, ~] = satMoonsConstants(4); mu_TBs(i) = mu; % Rhea
        case 606, [~, mu, ~, ~] = satMoonsConstants(5); mu_TBs(i) = mu; % Titan
        case 607, mu_TBs(i) = 0.374;   % Hyperion
        case 608, mu_TBs(i) = 120.4;   % Iapetus
        otherwise
            error('lc_vs_nbody_sma:UnknownNaifID', ...
                  'The NAIF ID %d is not recognized for mu calculation.', id);
    end
end
      
%--- Setup Ephemeris Handle and Propagate ---
ephem_handle = @(t_sec) get_all_moon_positions(t_sec, t0_mjd, actualBodyNaifIDs, spiceParam);

% --- Pass the special perturbation IDs to the propagator ---
[~, nbody_state] = propagateNBodyODE3(r0_sc', v0_sc', time_vector, ...
                                   mu_central_body, mu_TBs, ephem_handle, specialPerturbationIDs,[]);

nbody_r_norms = vecnorm(nbody_state(:, 1:3), 2, 2);
[ra_nbody, ~] = max(nbody_r_norms);
[rp_nbody, ~] = min(nbody_r_norms);

a_nbody = (ra_nbody + rp_nbody) / 2;

end