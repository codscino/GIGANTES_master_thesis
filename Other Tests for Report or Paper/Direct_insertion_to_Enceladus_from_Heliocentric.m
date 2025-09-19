clear all; close all; clc; format short g; warning off;
%% TESTCASE 2. PARTIAL COT Equinox Case
%--------------------------------------------------------------------------
% Here we search for a sequence of 10 consecutive fly-bys in Enceladus.
% This is similar to those proposed in the L4 report.
%% 1. Loading Constants & parameters
% --------------------------------------------------------------------------
% Constants & Parameters
RAD=pi/180;
Days2Sec=3600*24;
%--------------------------------------------------------------------------
% Define Central Body & Moon of interest
parsPLOT.INPUTS.idCentral  = 6; % (5 = Jupiter; 6 = Saturn; 7 = Uranus)
parsPLOT.INPUTS.idMoon     = 1; % (1 = Io; 2 = Europa; 3 = Ganymede; 4 = Callisto)
%(1 = Enceladus; 2 = Thetys; 3 = Dione; 4 = Rhea; 5 = Titan)
%(1 = Miranda; 2 = Ariel; 3 = Umbriel; 4 = Titania; 5 = Oberon)

[parsPLOT.Planet.mu, parsPLOT.Planet.EquRad, parsPLOT.Planet.OrbRad, parsPLOT.Planet.hmin] = planetConstants(parsPLOT.INPUTS.idCentral); %[km3/s2],[km],[km] & [km]

%Retrieve Desired Moon Parameters
if parsPLOT.INPUTS.idCentral == 3
    parsPLOT.Moon.OrbRad = 384748; parsPLOT.Moon.mu  = getAstroConstants('Moon','mu'); %[km],[km3/s2]
    parsPLOT.Moon.EquRad = getAstroConstants('Moon','Radius'); parsPLOT.Moon.hmin = 50;  %[km], [km]
elseif parsPLOT.INPUTS.idCentral == 5
    [parsPLOT.Moon.OrbRad, parsPLOT.Moon.mu, parsPLOT.Moon.EquRad, parsPLOT.Moon.hmin] = jupMoonsConstants(parsPLOT.INPUTS.idMoon); %[km],[km3/s2],[km] & [km]
elseif parsPLOT.INPUTS.idCentral == 6
    [parsPLOT.Moon.OrbRad, parsPLOT.Moon.mu, parsPLOT.Moon.EquRad, parsPLOT.Moon.hmin] = satMoonsConstants(parsPLOT.INPUTS.idMoon); %[km],[km3/s2],[km] & [km]
else
    [parsPLOT.Moon.OrbRad, parsPLOT.Moon.mu, parsPLOT.Moon.EquRad, parsPLOT.Moon.hmin] = uranusMoonsConstants(pars.INPUTS.idMoon); %[km],[km3/s2],[km] & [km]
end

for i = 1:size(parsPLOT.INPUTS.idMoon,2)
    parsPLOT.Moon.Vel(i)    = sqrt(parsPLOT.Planet.mu/parsPLOT.Moon.OrbRad(i));           %[km/s] Moon Orbital velocity
    parsPLOT.Moon.Period(i) = 2*pi*sqrt(parsPLOT.Moon.OrbRad(i)^3/parsPLOT.Planet.mu);    %[s] Moon orbital period
    parsPLOT.Moon.HillSph(i) = parsPLOT.Moon.OrbRad(i)*( parsPLOT.Moon.mu(i)/(3*(parsPLOT.Moon.mu(i) + parsPLOT.Planet.mu)))^(1/3);     %[km] Moon Hill's Sphere
end

%--------------------------------------------------------------------------
%Load some necessary parameters
%Define parameters regarding the flyby
parsPLOT.INPUTS.Flyby.min_h    = 25;   %[km] Minimum flyby altitude
parsPLOT.GroundTr.npoints      = 30e3; % N° of time steps for ground-track propagation
parsPLOT.GroundTr.t_prop       = 10;    %[minutes] Time of flyby hyperbola propagation
parsPLOT.INPUTS.Flyby.hMapping = 1500;  %[km] Max altitude to consider mapping

%Define starting epoch
parsPLOT.INPUTS.epoch0 = 0;
%Load sectirs for groundtrack plots
load('sectorObj.mat');
load('sectorObj_.mat');
parsPLOT.sectorObj  = sectorObj;
parsPLOT.sectorObj_ = sectorObj_;
parsPLOT.rect       = defineRectangleMapping(); % --> save the rectangle mapping
% --------------------------------------------------------------------------

VInf=norm([6.454 -1.118 0]);


muSaturn=parsPLOT.Planet.mu;

muEnceladus=parsPLOT.Moon.mu;
Sma_Enceladus=parsPLOT.Moon.OrbRad; 
RadiusEnceladus=parsPLOT.Moon.EquRad;


% Velocity of the Spacecraft at Arrival at Enceladus SOI.

Vp_Saturn=sqrt(2*(muSaturn/Sma_Enceladus+VInf^2/2));

vc_Enceladus=sqrt(muSaturn/Sma_Enceladus);

Vinf_min=Vp_Saturn-vc_Enceladus


% velocity of Spacecraft at 200 km of Enceladus Surface

Vp_Enceladus=sqrt(2*(muEnceladus/(RadiusEnceladus+200)+Vinf_min^2/2));


vc=sqrt(muEnceladus/(RadiusEnceladus+200))

DV_capture=Vp_Enceladus-vc

