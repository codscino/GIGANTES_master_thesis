clc;
clear all;
close all;

%% 1. Define Inputs
% ------------------
% Flyby parameters
pars.INPUTS.idCentral = 6;      % Saturn
pars.INPUTS.idMoon    = 1;      % Enceladus flyby
pars.INPUTS.V_inf     = 4.0;    % km/s
pars.INPUTS.alpha     = 0.15;   % rad
pars.INPUTS.k         = 0;

% SPICE and Time parameters
kernels = {'sat441.bsp', 'naif0012.tls', 'pck00010.tpc', 'de440.bsp'};
% kernels = {'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
loadSpiceKernels(kernels)
spiceParam.frame    = 'J2000';
spiceParam.abcorr   = 'NONE';
spiceParam.observer = '699'; % Saturn NAIF ID
t0_mjd = date2mjd2000([2025, 1, 1, 0, 0, 0]);

% N-Body Perturbers
% perturbingBodyNaifIDs = [602, 603, 605, 604, 601, 606, 10]; % Enceladus, Tethys, Rhea, Dione, Mimas, Titan, Sun
perturbingBodyNaifIDs = 606; %Titan


%% 2. Call the Function
[apoapsis_lc, apoapsis_kep, apoapsis_nbody,~] = ...
    calculateApoapsisComparison(pars, spiceParam, t0_mjd, perturbingBodyNaifIDs, 'backward', true);

[apoapsis_lc2, apoapsis_kep2, apoapsis_nbody2,~] = ...
    calculateApoapsisComparison(pars, spiceParam, t0_mjd, perturbingBodyNaifIDs, 'backward', true);


% cspice_kclear(); 

%% 3. Results
fprintf('\n--- Comparison of Apoapsis Results ---\n');
fprintf('Linked Conics: %.2f km\n', apoapsis_lc);
fprintf('Keplerian:     %.2f km\n', apoapsis_kep);
fprintf('N-Body:        %.2f km\n\n', apoapsis_nbody);

err_keplerian = apoapsis_nbody - apoapsis_kep;
err_linked_conics = apoapsis_nbody - apoapsis_lc;

fprintf('N-Body vs. Keplerian Error: %.2f km\n', err_keplerian);
fprintf('(This error is due to ignoring all third-body perturbations during propagation)\n\n');

fprintf('N-Body vs. Linked Conics Error: %.2f km\n', err_linked_conics);
fprintf('(This error is due to the idealized, analytical flyby calculation)\n');


%%
aa = linspace(0, 2*pi,36);
n = numel(aa);
f_in = zeros(n,1);

a_lc = zeros(n,1);
a_kep = zeros(n,1);
a_nbody = zeros(n,1);

err_lc = zeros(n,1);
err_kep = zeros(n,1);

for i = 1:n
    pars.INPUTS.alpha = aa(i);

    [a_lc(i), a_kep(i), a_nbody(i), f_in(i)] = calculateApoapsisComparison( ...
        pars, spiceParam, t0_mjd, perturbingBodyNaifIDs, 'forward', true);
    
    err_lc(i) = abs(a_nbody(i)-a_lc(i));
    err_kep(i) = abs(a_nbody(i)-a_kep(i));


end

% Plot everything together
figure;
plot(aa, f_in,   'b-',  'LineWidth', 1.5); hold on;
plot(aa, err_lc, 'r--', 'LineWidth', 1.5);
plot(aa, err_kep,'g-.','LineWidth', 1.5);
hold off;

xlabel('\alpha [rad]');
ylabel('Value');
title('Incoming Angle vs. f_{in} and Errors');
legend('f_{in}','err_{lc}','err_{kep}', 'Location','best');
grid on;