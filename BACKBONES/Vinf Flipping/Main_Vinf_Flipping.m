%% INITIALIZATION %%
clear all; close all; clc; format long g;
% addpath(genpath([pwd '/functions']));

% Define Central Body & Moon of interest
pars.INPUTS.idCentral  = 5; % (5 = Jupiter; 6 = Saturn; 7 = Uranus)
pars.INPUTS.idMoon     = [2]; % (1 = Io; 2 = Europa; 3 = Ganymede; 4 = Callisto)
                            % (1 = Enceladus; 2 = Thetys; 3 = Dione; 4 = Rhea; 5 = Titan)
                            % (1 = Miranda; 2 = Ariel; 3 = Umbriel; 4 = Titania; 5 = Oberon)

% Define Flyby Hyperbolic Excess Velocities
% pars.INPUTS.V_inf = 4.0;     %[km/s] for Enceladus 

% Define the crank angle to consider
pars.INPUTS.crank = [0, pi];

% Define the resonances to consider
pars.INPUTS.Resonances = [3 1];

% Define parameters regarding the flyby
pars.INPUTS.Flyby.min_h    = 25;   %[km] Minimum flyby altitude
pars.GroundTr.npoints      = 15e3; % N° of time steps for ground-track propagation
pars.GroundTr.t_prop       = 10;    %[minutes] Time of flyby hyperbola propagation 
pars.INPUTS.Flyby.hMapping = 300;  %[km] Max altitude to consider mapping 
pars.INPUTS.maxDV_Defect  = 0; %[km/s] Maximum DV defect allowed for a flyby

% Define starting epoch
pars.INPUTS.epoch0 = 0;

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

% Define other parameters
pars.AU     = 1.49597870691*10^(8);   %[km]
pars.g0     = 9.80665;                %[m/s^2] 
pars.Day    = 86400;                  %[s]
pars.Year   = 365.25;                 %[days]
pars.JD     = 2400000.5;

% % Determine maximum bending due to flyby
rp_flyby  = pars.INPUTS.Flyby.min_h + pars.Moon.EquRad;           %[km]
% V_infs = pars.INPUTS.V_inf;
% for i = 1:length(V_infs)
%     e_fly     = 1 + ((rp_flyby*V_infs(i)^2)/pars.Moon.mu); %[-] 
%     delta_max(i) = 2*asin(1/e_fly);                               %[rad]
% end
% pars.delta_max = delta_max;
pars.rp_flyby = rp_flyby;

% Load sectors for groundtrack plots
load('sectorObj.mat');
load('sectorObj_.mat');
pars.sectorObj  = sectorObj;
pars.sectorObj_ = sectorObj_;
pars.rect       = defineRectangleMapping(); % --> save the rectangle mapping

%% COMPUTE Vinf FLIPPING SOLUTION FOR THE SPECIFIC RESONANCE %%
% Define n° of flybys (for V-inf flipping this is 1)
pars.N_flybys = 1;

% Define Initial guess [Vinf (km/s); Alfa (rad)]
x0 = [3.2, 0]; 

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

    % Define the Vinf Flipping Outbound to Inbound Transfer
    VinfFlip_TransfersDatabase = [Vinf_Flip_Solutions(i).N,  Vinf_Flip_Solutions(i).M,  Vinf_Flip_Solutions(i).V_inf, ...
         Vinf_Flip_Solutions(i).Alfa, 0, Reso(1)*pars.Moon.Period, Vinf_Flip_Solutions(i).V_inf,  Vinf_Flip_Solutions(i).Alfa, pi, 810];  %[810 = 81 outbound/inbound + 0 because its a resonance]

    % Define the Vinf Flipping Inbound to Outbound Transfer
    VinfFlip_TransfersDatabase = [VinfFlip_TransfersDatabase; Vinf_Flip_Solutions(i).N,  Vinf_Flip_Solutions(i).M,  Vinf_Flip_Solutions(i).V_inf, ...
         Vinf_Flip_Solutions(i).Alfa, pi, Reso(1)*pars.Moon.Period, Vinf_Flip_Solutions(i).V_inf,  Vinf_Flip_Solutions(i).Alfa, 0, 180]; %[180 = 18 inbound/outbound + 0 because its a resonance]

end

% Check if d_req is d_max
d_req = acos(cos(Vinf_Flip_Solutions.Alfa)*cos(Vinf_Flip_Solutions.Alfa) + sin(Vinf_Flip_Solutions.Alfa)*sin(Vinf_Flip_Solutions.Alfa)*cos(pi)); 
d_max = 2*asin(1/(1 + ((rp_flyby*Vinf_Flip_Solutions.V_inf^2)/pars.Moon.mu)));

% Define the Vinf velocity
pars.INPUTS.V_inf =  Vinf_Flip_Solutions.V_inf;
pars.delta_max = d_max;

%% COMPUTE DATABASE OF PSEUDO-RESONANT TRANSFERS %%
% Define which 81 pseudo-resonance to use
pars.INPUTS.remove81 = 1;  %[0] = add 1 rev to the outbound-inbound short case
                           %[1] = do not add one rev to the outbound-inbound short case

PseudoReso_TransfersDatabase = [];
for i = 1:size(Vinf_Flip_Solutions.V_inf , 2)
    [PseudoReso_Transfer, PseudoResoStruct] = PseudoResoTransfers_Generation(Vinf_Flip_Solutions(i).V_inf, pars);
    PseudoReso_TransfersDatabase = [PseudoReso_TransfersDatabase; PseudoReso_Transfer];
end

%% CONSTRUCT PETAL ROTATION USING V-INF FLIPPING %%
% Setup inputs for Vinf Flipping sequences (name PetalRot left because
% functions from Petal Rotation are being re-used)
pars.PetalRot.Starting_Reso  = pars.INPUTS.Resonances; % Resonance Orbit at Start of Petal Rotation
pars.PetalRot.Starting_Crank = [0];    % [deg] Crank angle at start of Vinf Flip sequence (0 = Outbound / 180° = Inbound )
pars.PetalRot.Starting_Vinf  = Vinf_Flip_Solutions.V_inf;     % [km/s] Vinf at Start of Vinf Flip sequence
pars.PetalRot.Num_Flybys     = 36;     % Maximum n° of flybys that the Petal Rotation can use
pars.PetalRot.Beam_Width     = 1000;  % N° of solutions to store at each Petal Rotation expansion step

% Retrieve Vinf Flipping Option
if pars.PetalRot.Starting_Crank == 0
    VinfFlip_TransfersDatabase = VinfFlip_TransfersDatabase(1,:);  
elseif pars.PetalRot.Starting_Crank == 180
    VinfFlip_TransfersDatabase = VinfFlip_TransfersDatabase(2,:);
end

[Vinf_Flip_Seqs, Vinf_Flip_FlybysData] = VinfFlipRotation_BuildUp(VinfFlip_TransfersDatabase, PseudoReso_TransfersDatabase, pars);

%% ANALYZE RESULTS %%
% Select the index of a solution to analyze
idx = 1;
VinfFlip_Sequence = cell2mat(Vinf_Flip_Seqs(idx));
VinfFlip_FlybyData = cell2mat(Vinf_Flip_FlybysData(idx));

% Plot the Petal Rotation Planet-Centric Trajectory
[Fig] = Plot_PetalRotation_Trajectory(VinfFlip_Sequence, VinfFlip_FlybyData, pars);

% Plot the Petal Rotation Ground-Track
for i = 2:size(VinfFlip_Sequence, 1)
    node_in = VinfFlip_Sequence(i-1,7:9);
    node_out = VinfFlip_Sequence(i, 3:5);
    Flyby_Data(i) = Flyby_BuildUp(node_in, node_out, pars);
end

Flyby_Data(1) = [];
colors = cool(length(Flyby_Data));
fig1 = figure( 'Color', [1 1 1] );
hold on;
plotTextureLatLong(pars.INPUTS.idMoon , pars.INPUTS.idCentral , 1);
plotSquares(pars, 1);  
axis normal;
for i = 1:size(Flyby_Data, 2)
    Plot_Flyby_GT(Flyby_Data(i), colors(i,:));
end

