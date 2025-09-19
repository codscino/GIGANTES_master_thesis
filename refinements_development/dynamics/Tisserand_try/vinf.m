
clc;
clear all;
close all;

%% ========================================================================
%  1. DEFINE INPUTS & INITIALIZE
%  ========================================================================

soi_multiplier = 64; 

pars.GroundTr.t_prop  = 50*2000;
pars.GroundTr.npoints = 30e3;

pars.INPUTS.perturbingBodyNaifIDs = [-2, -3, -4, 606, 605, 10, 603, 604, 601, 607, 608, 602, 5];

pars.INPUTS.idCentral = 6;
pars.INPUTS.idMoon    = 1;
pars.INPUTS.V_inf     = 4.0;
pars.INPUTS.Flyby.min_h = 25;

kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels);

% --- Load Gravitational Parameters ---
[pars.Planet.mu, ~, ~, ~] = planetConstants(pars.INPUTS.idCentral);
[pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = satMoonsConstants(pars.INPUTS.idMoon);
pars.Moon.HillSph = pars.Moon.OrbRad * (pars.Moon.mu / (3 * pars.Planet.mu))^(1/3);

% Dynamically build mu_TBs list
actualBodyNaifIDs = pars.INPUTS.perturbingBodyNaifIDs(pars.INPUTS.perturbingBodyNaifIDs >= 0);
mu_TBs = zeros(1, length(actualBodyNaifIDs));
for i = 1:length(actualBodyNaifIDs)
    id = actualBodyNaifIDs(i);
    switch id
        case 10,  mu_TBs(i) = getAstroConstants('Sun', 'Mu');
        case 5, mu_TBs(i) = getAstroConstants('Jupiter', 'Mu');
        case 601, [~, mu, ~, ~] = satMoonsConstants(0); mu_TBs(i) = mu;
        case 602, [~, mu, ~, ~] = satMoonsConstants(1); mu_TBs(i) = mu;
        case 603, [~, mu, ~, ~] = satMoonsConstants(2); mu_TBs(i) = mu;
        case 604, [~, mu, ~, ~] = satMoonsConstants(3); mu_TBs(i) = mu;
        case 605, [~, mu, ~, ~] = satMoonsConstants(4); mu_TBs(i) = mu;
        case 606, [~, mu, ~, ~] = satMoonsConstants(5); mu_TBs(i) = mu;
        case 607, mu_TBs(i) = 0.374;
        case 608, mu_TBs(i) = 120.4;
    end
end
pars.INPUTS.mu_TBs = mu_TBs;


% --- Set a specific flyby epoch ---
pars.INPUTS.epoch0 = date2mjd2000([2035, 1, 10, 0, 0, 0]);

pars.INPUTS.NAIFCentral = 699;
pars.INPUTS.NAIFMoon = 602;


% propagate to 64 SOI (stable sma 3 body)
m_enceladus = 1.08e20;
m_saturn = getAstroConstants('Saturn', 'Mass');
a_enceladus = 238020;
r_soi_enceladus = a_enceladus * (m_enceladus / m_saturn)^(2/5);

pars.INPUTS.maxPropagationDistance = r_soi_enceladus*soi_multiplier;

%% ========================================================================
%  2. LOOP THROUGH DIFFERENT FLYBY GEOMETRIES
%  ========================================================================

% --- Define the range for the parameter to vary ---
v_inf_range = linspace(3, 5, 20);
num_sims = length(v_inf_range);

% --- Pre-allocate results arrays ---
lc_out_norms = zeros(num_sims, 1);
nb_out_norms = zeros(num_sims, 1);
deviation_at_exit = zeros(num_sims, 1);


tic;

for i = 1:num_sims
    
    pars.INPUTS.V_inf = v_inf_range(i);

    nodein = [pars.INPUTS.V_inf, 0.15, 0];
    nodeout = [pars.INPUTS.V_inf, 0.15, deg2rad(1)];

    % --- Call the main function ---
    [~, ~, ~, LC_out, NB_out] = propagateAndCompare(nodein, nodeout, pars);

    % --- Store the results ---
    lc_out_norms(i) = norm(LC_out(1:3));
    nb_out_norms(i) = norm(NB_out(1:3));
    deviation_at_exit(i) = norm(LC_out(1:3) - NB_out(1:3));
    
end

toc; 

%% ========================================================================
%  3. PLOT THE RESULTS
%  ========================================================================

