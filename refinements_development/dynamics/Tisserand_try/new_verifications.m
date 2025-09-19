clc;
clear all;
close all;

%% ========================================================================
%  1. DEFINE INPUTS & INITIALIZE
%  ========================================================================

% soi_multiplier = 136.6; 
soi_multiplier = 64; 

pars.GroundTr.t_prop  = 50*2000;
pars.GroundTr.npoints = 30e3;

pars.INPUTS.perturbingBodyNaifIDs = [-2,10,602];
% pars.INPUTS.perturbingBodyNaifIDs = [];


pars.INPUTS.idCentral = 6;
pars.INPUTS.idMoon    = 1;
pars.INPUTS.Flyby.min_h = 25;

kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels);

% --- Load Gravitational Parameters ---
[pars.Planet.mu, ~, ~, ~] = planetConstants(pars.INPUTS.idCentral);
[pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = satMoonsConstants(pars.INPUTS.idMoon);
pars.Moon.HillSph = pars.Moon.OrbRad * (pars.Moon.mu / (3 * pars.Planet.mu))^(1/3);

% Dynamically build mu_TBs list
actualBodyNaifIDs_init = pars.INPUTS.perturbingBodyNaifIDs(pars.INPUTS.perturbingBodyNaifIDs >= 0);
mu_TBs = zeros(1, length(actualBodyNaifIDs_init));
for i = 1:length(actualBodyNaifIDs_init)
    id = actualBodyNaifIDs_init(i);
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

pars.INPUTS.NAIFCentral = 699;
pars.INPUTS.NAIFMoon = 602;

% propagate to 64 SOI
m_enceladus = 1.08e20;
m_saturn = getAstroConstants('Saturn', 'Mass');
a_enceladus = 238020;
r_soi_enceladus = a_enceladus * (m_enceladus / m_saturn)^(2/5);
pars.INPUTS.maxPropagationDistance = r_soi_enceladus * soi_multiplier;

% --- Set baseline constant parameters for loops ---
BASELINE_VINF = 4;

pars.INPUTS.V_inf = BASELINE_VINF;

%% ========================================================================
%  2A. EPOCH LOOP (MUST RUN FIRST TO FIND OPTIMAL EPOCH)
%  ========================================================================

% Use a non-zero deflection for a realistic scenario
nodein_base = [BASELINE_VINF, 0, 0];
% nodeout_base = [BASELINE_VINF, 0.15, deg2rad(1)]; 
nodeout_base = [BASELINE_VINF, 0.15, deg2rad(1)]; 

pars.INPUTS.epoch0 = date2mjd2000([2040, 1, 1, 0, 0, 0]);

%% ====================================================================
%  START: Dismantled propagateAndCompare (Sections 1-3)
%  ====================================================================

%% Section 1: INITIALIZE PARAMETERS & SPICE
%  --------------------------------------------------------------------
spiceParam.frame    = 'J2000';
spiceParam.abcorr   = 'NONE';
spiceParam.observer = '699'; % Saturn Center NAIF ID

mu_central_body = pars.Planet.mu;
mu_enceladus = pars.Moon.mu;

% Perturbing bodies setup from the 'pars' struct
perturbingBodyNaifIDs = pars.INPUTS.perturbingBodyNaifIDs;
actualBodyNaifIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs >= 0);
specialPerturbationIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs < 0);
mu_TBs_local = pars.INPUTS.mu_TBs;



%% 3
[pLC, pNB] = calculate_pericenter_states(nodein_base, nodeout_base, pars);


%% 4
norm(pLC(1:3)-pNB(1:3))