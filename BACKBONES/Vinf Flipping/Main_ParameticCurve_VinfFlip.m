%% INITIALIZATION %%
clear all; close all; clc; format long g;
% addpath(genpath([pwd '/functions']));

% Define Central Body & Moon of interest
pars.INPUTS.idCentral  = 5; % (5 = Jupiter; 6 = Saturn; 7 = Uranus)
pars.INPUTS.idMoon     = [2]; % (1 = Io; 2 = Europa; 3 = Ganymede; 4 = Callisto)
                            % (1 = Enceladus; 2 = Thetys; 3 = Dione; 4 = Rhea; 5 = Titan)
                            % (1 = Miranda; 2 = Ariel; 3 = Umbriel; 4 = Titania; 5 = Oberon)

% Define Flyby Hyperbolic Excess Velocities
pars.INPUTS.V_inf = 4;     %[km/s] for Enceladus 

% Define the crank angle to consider
pars.INPUTS.crank = [0, pi];

% Define the resonances to consider
pars.INPUTS.Resonances = [5 1; 11 2; 6 1; 13 2; 7 1; 15 2; 8 1; 17 2];

% Define parameters regarding the flyby
pars.INPUTS.Flyby.min_h    = 25;   %[km] Minimum flyby altitude
pars.GroundTr.npoints      = 15e3; % N° of time steps for ground-track propagation
pars.GroundTr.t_prop       = 10;    %[minutes] Time of flyby hyperbola propagation 
pars.INPUTS.Flyby.hMapping = 300;  %[km] Max altitude to consider mapping 
pars.INPUTS.maxDV_Defect  = 0.1; %[km/s] Maximum DV defect allowed for a flyby

% Retrieve Central Body (Planet) Parameters
[pars.Planet.mu, pars.Planet.EquRad, pars.Planet.OrbRad, pars.Planet.hmin] = planetConstants(pars.INPUTS.idCentral); %[km3/s2],[km],[km] & [km]

% Retrieve Desired Moon Parameters 
if pars.INPUTS.idCentral == 3
    pars.Moon.OrbRad = 384748; pars.Moon.mu  = getAstroConstants('Moon','mu'); %[km],[km3/s2]
    pars.Moon.EquRad = getAstroConstants('Moon','Radius'); pars.Moon.hmin = 50;  %[km], [km]
elseif pars.INPUTS.idCentral == 5
    [pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = jupMoonsConstants(pars.INPUTS.idMoon); %[km],[km3/s2],[km] & [km]
elseif pars.INPUTS.idCentral == 6
    [pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = satMoonsConstants(pars.INPUTS.idMoon); %[km],[km3/s2],[km] & [km]
else
    [pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = uranusMoonsConstants(pars.INPUTS.idMoon); %[km],[km3/s2],[km] & [km]
end

for idx = 1:size(pars.INPUTS.idMoon,2)
    pars.Moon.Vel(idx)    = sqrt(pars.Planet.mu/pars.Moon.OrbRad(idx));           %[km/s] Moon Orbital velocity
    pars.Moon.Period(idx) = 2*pi*sqrt(pars.Moon.OrbRad(idx)^3/pars.Planet.mu);    %[s] Moon orbital period
    pars.Moon.HillSph(idx) = pars.Moon.OrbRad(idx)*( pars.Moon.mu(idx)/(3*(pars.Moon.mu(idx) + pars.Planet.mu)))^(1/3);     %[km] Moon Hill's Sphere
end

% Determine maximum bending due to flyby
rp_flyby  = pars.INPUTS.Flyby.min_h + pars.Moon.EquRad;     %[km]
V_infs = pars.INPUTS.V_inf;
for i = 1:length(V_infs)
    e_fly     = 1 + ((rp_flyby*V_infs(i)^2)/pars.Moon.mu);  %[-] 
    delta_max(i) = 2*asin(1/e_fly);                         %[rad]
end
pars.delta_max = delta_max;
pars.rp_flyby = rp_flyby;

%% COMPUTE Vinf FLIPPING FOR SPECIFIC RESONANCES %%
% Define n° of flybys (for V-inf flipping this is 1)
pars.N_flybys = 1;

% Define Initial guess [Vinf (km/s); Alfa (rad)]
x0 = [2, 0]; 

% Create the output vectors
dim = size(pars.INPUTS.idMoon,1);
Vinf_Flip_Solutions = struct('N', cell(1,size(dim,1)), ...
    'M', cell(1,size(dim,1)), ...
    'V_inf', cell(1,size(dim,1)),...
    'Alfa',cell(1,size(dim,1)));

for i = 1:size(pars.INPUTS.Resonances, 1)

    % Retrieve the resonance ratio
    Reso = [pars.INPUTS.Resonances(i,1), pars.INPUTS.Resonances(i,2)];

    % Solve the system of nonlinear equations for each resonance
    eqt_handle = @(x) roots_Vinf_Flipping(x, Reso, pars.N_flybys, pars);
    solutions(i,:) = fsolve(eqt_handle,x0);

    % Store the results
    Vinf_Flip_Solutions(i).N = Reso(1);
    Vinf_Flip_Solutions(i).M = Reso(2);
    Vinf_Flip_Solutions(i).V_inf = solutions(i,1);   %[km/s]
    Vinf_Flip_Solutions(i).Alfa = solutions(i,2);    %[rad]

end

%% COMPUTE Vinf FLIPPING PARAMETRIC CURVE %%
% Define Initial guess [Vinf (km/s); Alfa (rad)]
x0 = [2, 0]; 

% Define vector of Moon periods
Moon_Periods = linspace(1, 10, 1000);

for i = 1:size(Moon_Periods, 2)

    % Retrieve the resonance ratio
    Reso = [Moon_Periods(i), 1];

    % Solve the system of nonlinear equations for each resonance
    eqt_handle = @(x) roots_Vinf_Flipping(x, Reso, 1, pars);
    solutions(i,:) = fsolve(eqt_handle,x0);
end


VinfFlip_Curve = figure('Color', [1 1 1]);clf;set(VinfFlip_Curve,'defaulttextinterpreter','latex');
hold on; grid on;
plot(solutions(:,1), Moon_Periods, 'LineWidth', 2);
hold on;
xlabel('Hyperbolic Excess Velocity [km/s]', 'FontSize', 18);
ylabel('Spacecraft Orbital Period [Moon Periods]', 'FontSize', 18);

% %% COMPUTE FULL PARAMETRIC CURVES %%
% % Define Initial guess [Vinf (km/s); Alfa (rad); SC Orbital Velocity (km/s)]
% x0 = [2, 0, pars.Moon.Vel]; 
% 
% % Define vector of Moon periods
% Moon_Periods = linspace(1, 10, 1000);
% 
% for i = 1:size(Moon_Periods, 2)
%     
%     Tsc = pars.Moon.Period*Moon_Periods(i); %[seconds]
% 
%     % Solve the system of nonlinear equations for each resonance
%     eqt_handle = @(x) roots_Vinf_Flipping_CHECK(x, Tsc, 10, pars);
%     solut = fsolve(eqt_handle,x0);
% 
%     if isreal(solut)
%         A=1;
%     end
% end


