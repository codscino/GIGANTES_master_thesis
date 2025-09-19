clc;
clear all;
close all;
% script concerning flyby paramters and some check on them
%% ========================================================================
%  1. DEFINE INPUTS & INITIALIZE
%  ========================================================================

soi_multiplier = 64; 

pars.t_prop_hours  = 240; % hours
pars.GroundTr.npoints = 30e3;

pars.INPUTS.perturbingBodyNaifIDs = []; % no pertubers 

% Perturbing bodies setup from the 'pars' struct
perturbingBodyNaifIDs = pars.INPUTS.perturbingBodyNaifIDs;
specialPerturbationIDs = perturbingBodyNaifIDs(perturbingBodyNaifIDs < 0);

pars.INPUTS.idCentral = 6;
pars.INPUTS.idMoon    = 1;
pars.INPUTS.Flyby.min_h = 10;

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

pars.INPUTS.NAIFCentral = 699;
pars.INPUTS.NAIFMoon = 602;

% propagate to 64 SOI
m_enceladus = 1.08e20;
m_saturn = getAstroConstants('Saturn', 'Mass');
a_enceladus = 238020;
r_soi_enceladus = a_enceladus * (m_enceladus / m_saturn)^(2/5);
pars.INPUTS.maxPropagationDistance = r_soi_enceladus * soi_multiplier;

spiceParam.frame    = 'J2000';
spiceParam.abcorr   = 'NONE';
spiceParam.observer = '699'; % Saturn Center NAIF ID


%%% flyby paramters
pars.INPUTS.epoch0 = date2mjd2000([2040 1 1 12 0 0]);
pars.INPUTS.V_inf = 4;

%% ========================================================================
%  2. CALCULATIONS FOR PLOTS
%  ========================================================================

% --- Define the ranges ---
n_steps_2d = 200;
n_steps_3d = 50; % Fewer steps for the 3D plot to speed up calculation
alpha_out_range = linspace(0.0005, 0.0034, n_steps_2d);
k_out_range = linspace(0.0005, 0.0034, n_steps_2d);
alpha_out_range_3d = linspace(0.0005, 0.0034, n_steps_3d);
k_out_range_3d = linspace(0.00005, 0.0034, n_steps_3d);


% --- Pre-calculate constant values ---
[r_enceladus_at_flyby, v_enceladus_at_flyby] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, pars.INPUTS.epoch0, true, spiceParam);
Rm = buildRm(r_enceladus_at_flyby, v_enceladus_at_flyby);
e_fly     = 1 + (( (pars.INPUTS.Flyby.min_h + pars.Moon.EquRad) * pars.INPUTS.V_inf^2) / pars.Moon.mu);
delta_max = 2 * asin(1 / e_fly);

% --- Part A: Loop for subplot 1 (Vary alpha_out) ---
fprintf('Starting calculations for subplot 1...\n');
rp_values_alpha = NaN(1, n_steps_2d);
for i = 1:length(alpha_out_range)
    nodein = [pars.INPUTS.V_inf, 0, 0]; 
    nodeout = [pars.INPUTS.V_inf, alpha_out_range(i), 0.0005];
    [vvinfin, ~, ~, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodein(1), nodein(2), nodein(3), pars.INPUTS.epoch0, pars);
    [vvinfout, ~, ~, vvga] = vinfAlphaCrank_to_VinfCARTClaudio(nodeout(1), nodeout(2), nodeout(3), pars.INPUTS.epoch0, pars);
    [vvinfouBM, delta] = infVelBeforeDefect(vvinfin, vvinfout, delta_max, vvga);
    vvinfin_bf = (Rm * vvinfin')';
    Energy = 0.5 * norm(vvinfin_bf)^2;
    sma = -pars.Moon.mu / (2 * Energy);
    ecc = 1 / (sin(delta / 2));
    rp = sma * (1 - ecc);
    if isreal(rp) && rp >= 0, rp_values_alpha(i) = rp; end
end
fprintf('Finished.\n');

% --- Part B: Loop for subplot 2 (Vary k_out) ---
fprintf('Starting calculations for subplot 2...\n');
rp_values_k = NaN(1, n_steps_2d);
for i = 1:length(k_out_range)
    nodein = [pars.INPUTS.V_inf, 0, 0]; 
    nodeout = [pars.INPUTS.V_inf, 0.0005, k_out_range(i)];
    [vvinfin, ~, ~, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodein(1), nodein(2), nodein(3), pars.INPUTS.epoch0, pars);
    [vvinfout, ~, ~, vvga] = vinfAlphaCrank_to_VinfCARTClaudio(nodeout(1), nodeout(2), nodeout(3), pars.INPUTS.epoch0, pars);
    [vvinfouBM, delta] = infVelBeforeDefect(vvinfin, vvinfout, delta_max, vvga);
    vvinfin_bf = (Rm * vvinfin')';
    Energy = 0.5 * norm(vvinfin_bf)^2;
    sma = -pars.Moon.mu / (2 * Energy);
    ecc = 1 / (sin(delta / 2));
    rp = sma * (1 - ecc);
    if isreal(rp) && rp >= 0, rp_values_k(i) = rp; end
end
fprintf('Finished.\n');

% --- Part C: Loop for subplot 3 (Vary alpha_out and k_out) ---
fprintf('Starting calculations for subplot 3...\n');
rp_values_3d = NaN(length(alpha_out_range_3d), length(k_out_range_3d));
for i = 1:length(alpha_out_range_3d)
    for j = 1:length(k_out_range_3d)
        nodein = [pars.INPUTS.V_inf, 0, 0]; 
        nodeout = [pars.INPUTS.V_inf, alpha_out_range_3d(i), k_out_range_3d(j)];
        [vvinfin, ~, ~, ~] = vinfAlphaCrank_to_VinfCARTClaudio(nodein(1), nodein(2), nodein(3), pars.INPUTS.epoch0, pars);
        [vvinfout, ~, ~, vvga] = vinfAlphaCrank_to_VinfCARTClaudio(nodeout(1), nodeout(2), nodeout(3), pars.INPUTS.epoch0, pars);
        [vvinfouBM, delta] = infVelBeforeDefect(vvinfin, vvinfout, delta_max, vvga);
        vvinfin_bf = (Rm * vvinfin')';
        Energy = 0.5 * norm(vvinfin_bf)^2;
        sma = -pars.Moon.mu / (2 * Energy);
        ecc = 1 / (sin(delta / 2));
        rp = sma * (1 - ecc);
        if isreal(rp) && rp >= 0, rp_values_3d(i, j) = rp; end
    end
end
fprintf('Finished.\n');

%% ========================================================================
%  3. PLOT THE RESULTS
%  ========================================================================
figure;
% Make the figure window larger
set(gcf, 'Position', [100, 100, 1600, 500]);

% --- Subplot 1 ---
subplot(1, 3, 1);
plot(alpha_out_range, rp_values_alpha, 'LineWidth', 2);
title({'rp vs. Alpha Out', '(k_{out} = 0)'});
xlabel('Alpha Out (rad)');
ylabel('Pericenter Radius, rp (km)');
grid on;

% --- Subplot 2 ---
subplot(1, 3, 2);
plot(k_out_range, rp_values_k, 'LineWidth', 2, 'Color', 'r');
title({'rp vs. k Out', '(\alpha_{out} = 0)'});
xlabel('k Out (rad)');
ylabel('Pericenter Radius, rp (km)');
grid on;

% --- Subplot 3 ---
subplot(1, 3, 3);
surf(k_out_range_3d, alpha_out_range_3d, rp_values_3d', 'EdgeColor', 'none');
title('rp vs. Alpha Out and k Out');
xlabel('k Out (rad)');
ylabel('Alpha Out (rad)');
zlabel('Pericenter Radius, rp (km)');
colorbar;
view(3);