tic

%% STM with no deltaV
clc;
clear all;
close all;

%% ========================================================================
%  1. SETUP 
%  ========================================================================
soi_multiplier = 64; 

% Your existing parameter setup
pars.t_prop_hours = 45;
pars.GroundTr.npoints = 1000;
pars.INPUTS.perturbingBodyNaifIDs = [-2, 10, 602];
% pars.INPUTS.perturbingBodyNaifIDs = [-2,10,602]; % J2, Enceladus, Sun
% pars.INPUTS.perturbingBodyNaifIDs = [ -2, 606, 605, 10, 603, 604, 601, 607, 608, 5];
% pars.INPUTS.perturbingBodyNaifIDs = [602];
pars.INPUTS.idCentral = 6;
pars.INPUTS.idMoon = 1;
pars.INPUTS.Flyby.min_h = 10;

% Load kernels and constants (your existing code)
kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels);

[pars.Planet.mu, ~, ~, ~] = planetConstants(pars.INPUTS.idCentral);
[pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = satMoonsConstants(pars.INPUTS.idMoon);
pars.Moon.HillSph = pars.Moon.OrbRad * (pars.Moon.mu / (3 * pars.Planet.mu))^(1/3);

% Your other setup parameters
pars.INPUTS.NAIFCentral = 699;
pars.INPUTS.NAIFMoon = 602;

% propagate to 64 SOI (stable sma 3 body)
m_enceladus = 1.08e20;
m_saturn = getAstroConstants('Saturn', 'Mass');
a_enceladus = 238020;
r_soi_enceladus = a_enceladus * (m_enceladus / m_saturn)^(2/5);

pars.INPUTS.maxPropagationDistance = r_soi_enceladus*soi_multiplier;


pars.INPUTS.epoch0 = date2mjd2000([2030 1 1 0 0 0]);


pars.INPUTS.V_inf = 4;

% Flyby parameters
% nodein = [pars.INPUTS.V_inf, 0, 0];
% nodeout = [pars.INPUTS.V_inf, 0.15, deg2rad(1)];

% nodein =  [pars.INPUTS.V_inf, 1.5, 0];            % [km/s, rad, rad]
% nodeout = [pars.INPUTS.V_inf, 1.5, deg2rad(1)];   % [km/s, rad, rad]

nodein  = [4, deg2rad(45.6918), deg2rad(0)];
nodeout = [4, deg2rad(45.6918), deg2rad(1)];


%% ========================================================================
%  3. RUN STM OPTIMIZATION
%  ========================================================================

fprintf('\n========================================\n');
fprintf('STARTING STM OPTIMIZATION\n');
fprintf('========================================\n');

% Call the STM optimization function
[NB_in_corrected, B_achieved_final] = performSTMOptimization(pars, nodein, nodeout);

fprintf('\n========================================\n');
fprintf('STM OPTIMIZATION COMPLETE\n');
fprintf('========================================\n');

toc