clc;
clear all;
close all;

%% ========================================================================
%  1. DEFINE INPUTS & INITIALIZE
%  ========================================================================

soi_multiplier = 64; 

pars.GroundTr.t_prop  = 50*2000;
pars.GroundTr.npoints = 30e3;

% pars.INPUTS.perturbingBodyNaifIDs = []; 
pars.INPUTS.perturbingBodyNaifIDs = [-2, 602, 10];

pars.INPUTS.idCentral = 6;
pars.INPUTS.idMoon    = 1;
pars.INPUTS.Flyby.min_h = 25;

kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
% PARALLEL EDIT: Load kernels once for the main thread (pre-loops)
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

% PARALLEL EDIT: Start the parallel pool
if isempty(gcp('nocreate'))
    parpool(); % Starts the default parallel pool
end


%% ========================================================================
%  2A. EPOCH LOOP (MUST RUN FIRST TO FIND OPTIMAL EPOCH)
%  ========================================================================
disp('--- (1/6) Starting Epoch Loop to find minimum error epoch ---');
num_sims_epoch = 50;
date1 = date2mjd2000([2040, 1, 1, 0, 0, 0]);
date2 = date2mjd2000([2045, 1, 1, 0, 0, 0]);
epoch_range = linspace(date1, date2, num_sims_epoch);
pars.INPUTS.V_inf = BASELINE_VINF;

deltaT_values_epoch = zeros(num_sims_epoch, 1);
pos_diff_out_epoch = zeros(num_sims_epoch, 1);
pos_diff_pericenter_epoch = zeros(num_sims_epoch, 1);

nodein_base = [BASELINE_VINF, 0, 0];
nodeout_base = [BASELINE_VINF, deg2rad(0.01), deg2rad(0.01)]; 

tic;
% PARALLEL EDIT: Changed 'for' to 'parfor'
parfor i = 1:num_sims_epoch
    % PARALLEL EDIT: Load kernels for each worker
    loadSpiceKernels(kernels); 
    
    % Use a temporary 'pars' struct for this iteration to avoid communication overhead
    local_pars = pars;
    
    fprintf('Running Epoch sim %d of %d...\n', i, num_sims_epoch);
    local_pars.INPUTS.epoch0 = epoch_range(i);
    try
        [pLC, pNB, ~, LCout, NBout, Tin, Tout] = propagateAndCompare(nodein_base, nodeout_base, local_pars);
        deltaT_values_epoch(i) = abs(Tout - Tin);
        pos_diff_out_epoch(i) = norm(LCout(1:3) - NBout(1:3));
        pos_diff_pericenter_epoch(i) = norm(pLC(1:3) - pNB(1:3));
    catch ME
        fprintf('Error in Epoch sim %d: %s\n', i, ME.message);
        deltaT_values_epoch(i)=NaN; pos_diff_out_epoch(i)=NaN; pos_diff_pericenter_epoch(i)=NaN;
    end
end

[~, idx_min] = min(pos_diff_pericenter_epoch);
min_epoch = epoch_range(idx_min);
deltaT_at_min_epoch = deltaT_values_epoch(idx_min);
pos_diff_out_at_min_epoch = pos_diff_out_epoch(idx_min);

fprintf('Epoch loop finished in %.2f sec.\n', toc);
fprintf('Minimum pericenter error found at epoch MJD: %.2f\n\n', min_epoch);


%% ========================================================================
%  2B. V_INFINITY LOOP (RUNS AT THE OPTIMAL min_epoch)
%  ========================================================================
disp('--- (2/6) Starting V_infinity Loop (at min_epoch) ---');
num_sims_vinf = 50; 
v_inf_range = linspace(3.0, 6.0, num_sims_vinf);
pars.INPUTS.epoch0 = min_epoch;

deltaT_values_vinf = zeros(num_sims_vinf, 1);
pos_diff_out_vinf = zeros(num_sims_vinf, 1);
pos_diff_pericenter_vinf = zeros(num_sims_vinf, 1);

tic;
% PARALLEL EDIT: Changed 'for' to 'parfor'
parfor i = 1:num_sims_vinf
    % PARALLEL EDIT: Load kernels for each worker
    loadSpiceKernels(kernels);

    local_pars = pars;
    
    fprintf('Running V_inf sim %d of %d...\n', i, num_sims_vinf);
    local_pars.INPUTS.V_inf = v_inf_range(i);
    nodein = [v_inf_range(i), 0, 0];
    nodeout = [v_inf_range(i), deg2rad(0.01), deg2rad(0.01)]; 
    try
        [pLC, pNB, ~, LCout, NBout, Tin, Tout] = propagateAndCompare(nodein, nodeout, local_pars);
        deltaT_values_vinf(i) = abs(Tout - Tin);
        pos_diff_out_vinf(i) = norm(LCout(1:3) - NBout(1:3));
        pos_diff_pericenter_vinf(i) = norm(pLC(1:3) - pNB(1:3));
    catch ME
        fprintf('Error in V_inf sim %d: %s\n', i, ME.message);
        deltaT_values_vinf(i)=NaN; pos_diff_out_vinf(i)=NaN; pos_diff_pericenter_vinf(i)=NaN;
    end
end
fprintf('V_infinity loop finished in %.2f sec.\n\n', toc);


%% ========================================================================
%  2C. CALCULATE THE "ZERO" VALUE FOR 3D PLOT NORMALIZATION
%  ========================================================================
disp('--- (3/6) Calculating Zero-Point for Normalization ---');
pars.INPUTS.V_inf = BASELINE_VINF;
pars.INPUTS.epoch0 = min_epoch;
[pLC_z, pNB_z, ~, LCout_z, NBout_z, Tin_z, Tout_z] = propagateAndCompare(nodein_base, nodeout_base, pars);
deltaT_zero = abs(Tout_z - Tin_z);
pos_diff_out_zero = norm(LCout_z(1:3) - NBout_z(1:3));
pos_diff_pericenter_zero = norm(pLC_z(1:3) - pNB_z(1:3));
fprintf('Zero-point calculated.\n\n');


%% ========================================================================
%  2D. ALPHA_IN / ALPHA_OUT GRID (k=const)
%  ========================================================================
disp('--- (4/6) Starting alpha_in vs alpha_out Loop ---');
num_sims_grid = 20;
alpha_in_range = linspace(0, pi, num_sims_grid);
alpha_out_range = linspace(0, pi, num_sims_grid);

deltaT_alpha_grid = zeros(num_sims_grid, num_sims_grid);
pos_diff_out_alpha_grid = zeros(num_sims_grid, num_sims_grid);
pos_diff_pericenter_alpha_grid = zeros(num_sims_grid, num_sims_grid);
tic;
% PARALLEL EDIT: Changed outer 'for' to 'parfor'
parfor i = 1:num_sims_grid
    % PARALLEL EDIT: Load kernels for each worker
    loadSpiceKernels(kernels);

    local_pars = pars;
    
    for j = 1:num_sims_grid
        fprintf('Running alpha grid sim (%d, %d)...\n', i, j);
        nodein = [BASELINE_VINF, alpha_in_range(i), 0];
        nodeout = [BASELINE_VINF, alpha_out_range(j), deg2rad(0.01)];
        
        try
            [pLC, pNB, ~, LCout, NBout, Tin, Tout] = propagateAndCompare(nodein, nodeout, local_pars);
            deltaT_alpha_grid(i,j) = abs(Tout - Tin);
            pos_diff_out_alpha_grid(i,j) = norm(LCout(1:3) - NBout(1:3));
            pos_diff_pericenter_alpha_grid(i,j) = norm(pLC(1:3) - pNB(1:3));
        catch ME
            deltaT_alpha_grid(i,j)=NaN; pos_diff_out_alpha_grid(i,j)=NaN; pos_diff_pericenter_alpha_grid(i,j)=NaN;
        end
    end
end
fprintf('Alpha grid loop finished in %.2f sec.\n\n', toc);


%% ========================================================================
%  2E. K_IN / K_OUT GRID (alpha=const)
%  ========================================================================
disp('--- (5/6) Starting k_in vs k_out Loop ---');
k_in_range = linspace(0, pi, num_sims_grid);
k_out_range = linspace(0, pi, num_sims_grid);

deltaT_k_grid = zeros(num_sims_grid, num_sims_grid);
pos_diff_out_k_grid = zeros(num_sims_grid, num_sims_grid);
pos_diff_pericenter_k_grid = zeros(num_sims_grid, num_sims_grid);
tic;
% PARALLEL EDIT: Changed outer 'for' to 'parfor'
parfor i = 1:num_sims_grid
    % PARALLEL EDIT: Load kernels for each worker
    loadSpiceKernels(kernels);
    
    local_pars = pars;

    for j = 1:num_sims_grid
        fprintf('Running k grid sim (%d, %d)...\n', i, j);
        nodein = [BASELINE_VINF, 0, k_in_range(i)];
        nodeout = [BASELINE_VINF, deg2rad(0.01), k_out_range(j)];
        try
            [pLC, pNB, ~, LCout, NBout, Tin, Tout] = propagateAndCompare(nodein, nodeout, local_pars);
            deltaT_k_grid(i,j) = abs(Tout - Tin);
            pos_diff_out_k_grid(i,j) = norm(LCout(1:3) - NBout(1:3));
            pos_diff_pericenter_k_grid(i,j) = norm(pLC(1:3) - pNB(1:3));
        catch ME
            deltaT_k_grid(i,j)=NaN; pos_diff_out_k_grid(i,j)=NaN; pos_diff_pericenter_k_grid(i,j)=NaN; 
        end
    end
end
fprintf('K grid loop finished in %.2f sec.\n\n', toc);


%% ========================================================================
%  2F. EPOCH vs V_INF GRID (FOR 3D PLOTS)
%  ========================================================================
disp('--- (6/6) Starting Epoch vs V_inf Grid Loop ---');
num_sims_2d_grid = 25;
epoch_range_grid = linspace(date1, date2, num_sims_2d_grid);
v_inf_range_grid = linspace(3.0, 6.0, num_sims_2d_grid);

deltaT_epoch_vinf_grid = zeros(num_sims_2d_grid, num_sims_2d_grid);
pos_diff_out_epoch_vinf_grid = zeros(num_sims_2d_grid, num_sims_2d_grid);
pos_diff_pericenter_epoch_vinf_grid = zeros(num_sims_2d_grid, num_sims_2d_grid);

tic;
% PARALLEL EDIT: Changed outer 'for' to 'parfor'
parfor i = 1:num_sims_2d_grid % Loop for epoch
    % PARALLEL EDIT: Load kernels for each worker
    loadSpiceKernels(kernels);

    local_pars = pars;

    for j = 1:num_sims_2d_grid % Loop for V_inf
        fprintf('Running Epoch/V_inf grid sim (%d, %d)...\n', i, j);
        local_pars.INPUTS.epoch0 = epoch_range_grid(i);
        local_pars.INPUTS.V_inf = v_inf_range_grid(j);
        nodein = [v_inf_range_grid(j), 0, 0];
        nodeout = [v_inf_range_grid(j), deg2rad(0.01), deg2rad(0.01)];

        try
            [pLC, pNB, ~, LCout, NBout, Tin, Tout] = propagateAndCompare(nodein, nodeout, local_pars);
            deltaT_epoch_vinf_grid(i,j) = abs(Tout - Tin);
            pos_diff_out_epoch_vinf_grid(i,j) = norm(LCout(1:3) - NBout(1:3));
            pos_diff_pericenter_epoch_vinf_grid(i,j) = norm(pLC(1:3) - pNB(1:3));
        catch ME
            deltaT_epoch_vinf_grid(i,j)=NaN; 
            pos_diff_out_epoch_vinf_grid(i,j)=NaN; 
            pos_diff_pericenter_epoch_vinf_grid(i,j)=NaN;
        end
    end
end
fprintf('Epoch/V_inf grid loop finished in %.2f sec.\n\n', toc);


%% ========================================================================
%  3. PLOT THE RESULTS
%  ========================================================================
clc
disp('All simulations complete. Now plotting results...');
%% ========================================================================
%  3. PLOT THE RESULTS
%  ========================================================================

% --- FIGURE 1, 2, 3 (As before) ---
% Figure 1
figure('Name', 'Tisserand Parameter Change', 'Position', [50 50 1200 800]);
sgtitle('Sensitivity of Tisserand Parameter Change |\DeltaT|', 'FontSize', 16, 'FontWeight', 'bold');
subplot(2, 2, 1); plot(epoch_range, deltaT_values_epoch, 'r-', 'LineWidth', 1.5); grid on; hold on;
plot(min_epoch, deltaT_at_min_epoch, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
title('vs. Flyby Epoch (Raw)'); xlabel('Flyby Epoch (MJD2000)'); ylabel('|\DeltaT|'); legend('Epoch variation', 'Selected Min Error Epoch', 'Location', 'best');
subplot(2, 2, 2); plot(v_inf_range, deltaT_values_vinf - deltaT_at_min_epoch, 'b-', 'LineWidth', 1.5); grid on;
title('vs. V_{\infty} (Normalized by Epoch Contribution)'); xlabel('V_{\infty} (km/s)'); ylabel('Normalized |\DeltaT|');
[ALPHA_OUT_GRID, ALPHA_IN_GRID] = meshgrid(rad2deg(alpha_out_range), rad2deg(alpha_in_range));
subplot(2, 2, 3); surf(ALPHA_OUT_GRID, ALPHA_IN_GRID, deltaT_alpha_grid' - deltaT_zero); shading interp; colorbar;
title('vs. \alpha_{in}/\alpha_{out} (Normalized)'); xlabel('\alpha_{out} (deg)'); ylabel('\alpha_{in} (deg)'); zlabel('Normalized |\DeltaT|');
[K_OUT_GRID, K_IN_GRID] = meshgrid(rad2deg(k_out_range), rad2deg(k_in_range));
subplot(2, 2, 4); surf(K_OUT_GRID, K_IN_GRID, deltaT_k_grid' - deltaT_zero); shading interp; colorbar;
title('vs. k_{in}/k_{out} (Normalized)'); xlabel('k_{out} (deg)'); ylabel('k_{in} (deg)'); zlabel('Normalized |\DeltaT|');

% Figure 2
figure('Name', 'Position Difference at Exit', 'Position', [100 100 1200 800]);
sgtitle('Sensitivity of Position Difference at Exit Boundary', 'FontSize', 16, 'FontWeight', 'bold');
subplot(2, 2, 1); plot(epoch_range, pos_diff_out_epoch, 'r-', 'LineWidth', 1.5); grid on; hold on;
plot(min_epoch, pos_diff_out_at_min_epoch, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
title('vs. Flyby Epoch (Raw)'); xlabel('Flyby Epoch (MJD2000)'); ylabel('Position Difference (km)'); legend('Epoch variation', 'Selected Min Error Epoch', 'Location', 'best');
subplot(2, 2, 2); plot(v_inf_range, pos_diff_out_vinf - pos_diff_out_at_min_epoch, 'b-', 'LineWidth', 1.5); grid on;
title('vs. V_{\infty} (Normalized)'); xlabel('V_{\infty} (km/s)'); ylabel('Normalized Diff (km)');
subplot(2, 2, 3); surf(ALPHA_OUT_GRID, ALPHA_IN_GRID, pos_diff_out_alpha_grid' - pos_diff_out_zero); shading interp; colorbar;
title('vs. \alpha_{in}/\alpha_{out} (Normalized)'); xlabel('\alpha_{out} (deg)'); ylabel('\alpha_{in} (deg)'); zlabel('Normalized Diff (km)');
subplot(2, 2, 4); surf(K_OUT_GRID, K_IN_GRID, pos_diff_out_k_grid' - pos_diff_out_zero); shading interp; colorbar;
title('vs. k_{in}/k_{out} (Normalized)'); xlabel('k_{out} (deg)'); ylabel('k_{in} (deg)'); zlabel('Normalized Diff (km)');

% Figure 3
figure('Name', 'Position Difference at Pericenter', 'Position', [150 150 1200 800]);
sgtitle('Sensitivity of Position Difference at Pericenter', 'FontSize', 16, 'FontWeight', 'bold');
subplot(2, 2, 1); plot(epoch_range, pos_diff_pericenter_epoch, 'r-', 'LineWidth', 1.5); grid on; hold on;
plot(min_epoch, min(pos_diff_pericenter_epoch), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
title('vs. Flyby Epoch (Raw)'); xlabel('Flyby Epoch (MJD2000)'); ylabel('Position Difference (km)'); legend('Epoch variation', 'Selected Min Error Epoch', 'Location', 'best');
subplot(2, 2, 2); plot(v_inf_range, pos_diff_pericenter_vinf - min(pos_diff_pericenter_epoch), 'b-', 'LineWidth', 1.5); grid on;
title('vs. V_{\infty} (Normalized)'); xlabel('V_{\infty} (km/s)'); ylabel('Normalized Diff (km)');
subplot(2, 2, 3); surf(ALPHA_OUT_GRID, ALPHA_IN_GRID, pos_diff_pericenter_alpha_grid' - pos_diff_pericenter_zero); shading interp; colorbar;
title('vs. \alpha_{in}/\alpha_{out} (Normalized)'); xlabel('\alpha_{out} (deg)'); ylabel('\alpha_{in} (deg)'); zlabel('Normalized Diff (km)');
subplot(2, 2, 4); surf(K_OUT_GRID, K_IN_GRID, pos_diff_pericenter_k_grid' - pos_diff_pericenter_zero); shading interp; colorbar;
title('vs. k_{in}/k_{out} (Normalized)'); xlabel('k_{out} (deg)'); ylabel('k_{in} (deg)'); zlabel('Normalized Diff (km)');


% --- FIGURE 4: Combined Epoch and V_inf Sensitivity -- [NEW SECTION] --
figure('Name', 'Combined Epoch and V_inf Sensitivity', 'Position', [200 200 1200 800]);
sgtitle('Combined Sensitivity to Epoch and V_{\infty}', 'FontSize', 16, 'FontWeight', 'bold');
[VINF_GRID, EPOCH_GRID] = meshgrid(v_inf_range_grid, epoch_range_grid);

% --- Plot 1: Tisserand Change ---
subplot(2, 2, 1);
surf(VINF_GRID, EPOCH_GRID, deltaT_epoch_vinf_grid);
shading interp; colorbar;
title('|\DeltaT|');
xlabel('V_{\infty} (km/s)');
ylabel('Epoch (MJD2000)');
zlabel('|\DeltaT|');
view(3);

% --- Plot 2: Position Difference at Exit ---
subplot(2, 2, 2);
surf(VINF_GRID, EPOCH_GRID, pos_diff_out_epoch_vinf_grid);
shading interp; colorbar;
title('Position Difference at Exit Boundary');
xlabel('V_{\infty} (km/s)');
ylabel('Epoch (MJD2000)');
zlabel('Position Difference (km)');
view(3);

% --- Plot 3: Position Difference at Pericenter ---
subplot(2, 2, 3);
surf(VINF_GRID, EPOCH_GRID, pos_diff_pericenter_epoch_vinf_grid);
shading interp; colorbar;
title('Position Difference at Pericenter');
xlabel('V_{\infty} (km/s)');
ylabel('Epoch (MJD2000)');
zlabel('Position Difference (km)');
view(3);

% --- Turn off the last empty subplot ---
subplot(2,2,4);
axis off;