clear all; close all; clc;

%% campagnola params
camp_path = which('campagnola.mat');
load(camp_path);

% find the correct flyby
mask = strcmp(CAMPAGNOLA_TABLE.Flyby, '21E14'); % flyby 21E14
rowIdx = find(mask);

% extract epoch
date_str = CAMPAGNOLA_TABLE.DateET{rowIdx};
dt = datetime(date_str, 'InputFormat', 'dd-MMM.-yyyy HH:mm:ss');
time_arr = [ year(dt), month(dt), day(dt), hour(dt), minute(dt), second(dt) ];
epoch0 = date2mjd2000(time_arr);
pars.INPUTS.epoch0 = epoch0;

% nodein
vinfin = CAMPAGNOLA_TABLE.v_inf(rowIdx);
alfain = CAMPAGNOLA_TABLE.alpha_minus(rowIdx);
kin = CAMPAGNOLA_TABLE.kappa_minus(rowIdx);

pars.INPUTS.V_inf = vinfin;

% nodeout
alfaout = CAMPAGNOLA_TABLE.alpha_plus(rowIdx);
kout = CAMPAGNOLA_TABLE.kappa_plus(rowIdx);
% vinfout_x  = vinfin*sin(alfaout)*cos(kout);
% vinfout_y = vinfin*cos(alfaout);
% vinfout_z = -vinfin*sin(alfaout)*sin(kout);
% vinfout = sqrt(vinfout_x^2+vinfout_y^2+vinfout_z^2);
vinfout = vinfin;

% altitude
pars.INPUTS.Flyby.min_h = CAMPAGNOLA_TABLE.ALT(rowIdx);



%% DEFINE PLANET, MOON PARAMETERS & CONSTANTS %%
% Define Central Body & Moon of interest
pars.INPUTS.idCentral  = 5; % (5 = Jupiter; 6 = Saturn; 7 = Uranus)
pars.INPUTS.idMoon     = 2; % (1 = Io; 2 = Europa; 3 = Ganymede; 4 = Callisto)
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
% Define parameters regarding the flyby
pars.GroundTr.npoints      = 30e3; % N° of time steps for ground-track propagation
pars.GroundTr.t_prop       = 5;    %[minutes] Time of flyby hyperbola propagation  ???
pars.INPUTS.Flyby.hMapping = 300;  %[km] Max altitude to consider mapping ???

% Determine maximum bending due to flyby
rp_flyby  = pars.INPUTS.Flyby.min_h + pars.Moon.EquRad;           %[km
e_fly     = 1 + ((rp_flyby*pars.INPUTS.V_inf^2)/pars.Moon.mu);    %[-] 
delta_max = 2*asin(1/e_fly);                                      %[rad]
pars.delta_max = delta_max;


%% COMPUTE FLYBY CHARACTERISTICS %%
% Define incoming & outgoing nodes
nodein  = [vinfin, alfain, kin]; %[km/s, rad, rad]
nodeout = [vinfout, alfaout, kout];

% Compute flyby parameters
[Flyby] = Flyby_BuildUp(nodein, nodeout, pars);


%% compute differences: Campagnola vs GIGANTES
