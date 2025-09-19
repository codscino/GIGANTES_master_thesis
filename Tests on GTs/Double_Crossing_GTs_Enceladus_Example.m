clear all; close all; clc; warning off

% This script should generate and help visualizing a full fly-by with a
% given planet or moon. The idea is to generate both the ground track plot
% and the b-plane trajectory plot.

% obviously the script needs to input some change of node otherwise the
% fly-by is undefined. 

%% DEFINE PLANET, MOON PARAMETERS & CONSTANTS %%
% Define Central Body & Moon of interest
pars.INPUTS.idCentral  = 6; % (5 = Jupiter; 6 = Saturn; 7 = Uranus)
pars.INPUTS.idMoon     = [1]; % (1 = Io; 2 = Europa; 3 = Ganymede; 4 = Callisto)
                            % (1 = Enceladus; 2 = Thetys; 3 = Dione; 4 = Rhea; 5 = Titan)
                            % (1 = Miranda; 2 = Ariel; 3 = Umbriel; 4 = Titania; 5 = Oberon)

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

for i = 1:size(pars.INPUTS.idMoon,2)
    pars.Moon.Vel(i)    = sqrt(pars.Planet.mu/pars.Moon.OrbRad(i));           %[km/s] Moon Orbital velocity
    pars.Moon.Period(i) = 2*pi*sqrt(pars.Moon.OrbRad(i)^3/pars.Planet.mu);    %[s] Moon orbital period
    pars.Moon.HillSph(i) = pars.Moon.OrbRad(i)*( pars.Moon.mu(i)/(3*(pars.Moon.mu(i) + pars.Planet.mu)))^(1/3);     %[km] Moon Hill's Sphere
end

%% DEFINE FLYBY PARAMETERS & CONSTANTS %%
% Define Flyby Hyperbolic Excess Velocities
pars.INPUTS.V_inf = 4;     %[km/s] 

vinfin = pars.INPUTS.V_inf;


% Define parameters regarding the flyby
pars.INPUTS.Flyby.min_h    = 25;   %[km] Minimum flyby altitude
pars.GroundTr.npoints      = 30e3; % N° of time steps for ground-track propagation
pars.GroundTr.t_prop       = 30;    %[minutes] Time of flyby hyperbola propagation 
pars.INPUTS.Flyby.hMapping = 300;  %[km] Max altitude to consider mapping 

% Determine maximum bending due to flyby
rp_flyby  = pars.INPUTS.Flyby.min_h + pars.Moon.EquRad;           %[km
e_fly     = 1 + ((rp_flyby*pars.INPUTS.V_inf^2)/pars.Moon.mu);    %[-] 
delta_max = 2*asin(1/e_fly);                                      %[rad]
pars.delta_max = delta_max;

% Define starting epoch
pars.INPUTS.epoch0 = 0;

% Load sectirs for groundtrack plots
load('sectorObj.mat');
load('sectorObj_.mat');
pars.sectorObj  = sectorObj;
pars.sectorObj_ = sectorObj_;
pars.rect       = defineRectangleMapping(); % --> save the rectangle mapping

%% COMPUTE FLYBY CHARACTERISTICS %%

NODEIN=[4	0.349065850398866	1.57079632679490;...
        4	0.349065850398866	1.57079632679490;...
        4	0.349065850398866	4.71238898038469;...
        4	0.349065850398866	4.71238898038469];

NODEOUT=[4	0.349098647237303	1.56620484268886;...
         4	0.349098647237303	1.57538781090093;...
         4	0.349098647237303	4.70779749627865;...
         4	0.349098647237303	4.71698046449073];


%%

% Compute flyby parameters
for f=1:length(NODEIN(:,1))
Flyby(f) = Flyby_BuildUp_BC(NODEIN(f,:), NODEOUT(f,:), pars);
end


%% PLOT GROUDNTRACK %%
colors = cool(length(Flyby));
fig1 = figure( 'Color', [1 1 1] );
hold on;
plotTextureLatLong(pars.INPUTS.idMoon , pars.INPUTS.idCentral , 1);
plotSquares(pars, 1);  
axis normal;
for i = 1:size(Flyby, 2)
    Plot_Flyby_GT(Flyby(i), colors(i,:));
end
% 
% Flyby = Flyby_BuildUp_BC(NODEIN2, NODEOUT2, pars);
% for i = 1:size(Flyby, 2)
%     Plot_Flyby_GT(Flyby(i), colors(i,:));
% end

NODEIN_inDeg=[NODEIN(:,1) rad2deg(NODEIN(:,2)) rad2deg(NODEIN(:,3))]
NODEOUT_inDeg=[NODEOUT(:,1) rad2deg(NODEOUT(:,2)) rad2deg(NODEOUT(:,3))]