 clear all; close all; clc; format short g; warning off;
% script_Refine_Solutions.m
%--------------------------------------------------------------------------
load('Saturn_Enceladus_Titan_exploration_results_full4c.mat')
% load('Saturn_Enceladus_Titan_exploration_results_full4d.mat')
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
elseif pars.INPUTS.idCentral == 6
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
parsPLOT.GroundTr.t_prop       = 15;    %[minutes] Time of flyby hyperbola propagation
parsPLOT.INPUTS.Flyby.hMapping = 1500;  %[km] Max altitude to consider mapping

%Define starting epoch
parsPLOT.INPUTS.epoch0 = 0;
%Load sectirs for groundtrack plots
load('sectorObj.mat');
load('sectorObj_.mat');
parsPLOT.sectorObj  = sectorObj;
parsPLOT.sectorObj_ = sectorObj_;
parsPLOT.rect       = defineRectangleMapping(); % --> save the rectangle mapping
%%
indexMax_Score=find(Score_ScientificInteres==max(Score_ScientificInteres));

treeExploration_SubSet=treeExploration(indexMax_Score);
indexValue=[1:1:length(indexMax_Score)];
InfoSubSet=[indexValue' numEnceladusFB(indexMax_Score)' Delta_V(indexMax_Score)' Total_time_of_flight(indexMax_Score)' Delta_longitude(indexMax_Score)']


[~,indexSorted]=sort(InfoSubSet(:,3));

InfoSubSet=InfoSubSet(indexSorted,:);

    indexRelevantFlyBys=find(AFull_route(:,8)==1);
    length(indexRelevantFlyBys)
% ploting the maximum Diversity  tour
IndexCheck=8894  % This is the exemple to reanalize. 
AFull_route=treeExploration_SubSet{IndexCheck}.route;
totalDV=sum(treeExploration_SubSet{IndexCheck}.routeDVCosts)
routeAnatomy=treeExploration_SubSet{IndexCheck}.routeAnatomy;
plotfullSaturnTransfer(AFull_route,pars)

%% Refinement 
% this fuction must rereun the precomputed full route. Note the end node
% and starting node on the same planet may now have a small discontinutity
% due the fly-by itself. 
% since Lambert arc won't work, a MGA+DSM model is used, where one variable
% is optimized. This will be the fraction of time to the DSM, since
% everything else is already prefixed. 

Transfer2Refine.Full_route=AFull_route;

% First step it is to reconctruct the DV and full transfer. 
% --------------------------------------------------------------------------
% Begining of the function here
% --------------------------------------------------------------------------
% Constants & Parameters
RAD=pi/180;
Days2Sec=3600*24;
nEnceladus=sqrt(pars.Planet.mu/pars.Moon.OrbRad(1)^3);
nTitan=sqrt(pars.Planet.mu/pars.Moon.OrbRad(2)^3);
% --------------------------------------------------------------------------
% Epremeris Moons Set up
KepTitan=[pars.Moon.OrbRad(2) 0 0 0 0 0];
KepEnceladus=[pars.Moon.OrbRad(1) 0 0 0 0 0];
SVTitan_func = @(TrA)kep2car([KepTitan(1:5) TrA],pars.Planet.mu);
SVEnceladus_func = @(TrA)kep2car([KepEnceladus(1:5) TrA],pars.Planet.mu);
% --------------------------------------------------------------------------

Full_route=Transfer2Refine.Full_route;

% Step 1. Identify all NODES IN
indexEnceladus_NODEIN=find(Full_route(:,8)==1);
NODEIN_Matrix=Full_route(indexEnceladus_NODEIN,9:11);

% step 2. Identify all NODES OUT
indexEnceladus_NODEOUT=find(Full_route(:,1)==1);
NODEOUT_Matrix=Full_route(indexEnceladus_NODEOUT,2:4);

if length(indexEnceladus_NODEIN)>length(indexEnceladus_NODEOUT)
    NODEOUT_Matrix(end+1,:)=NODEIN_Matrix(end,:);
end

JumpsTrack=[0; diff(indexEnceladus_NODEIN)-1];

%--------------------------------------------------------------------------
% intermediate Plotting of Individual Fly-bys

indexNode=8;

run('Plot_singleFlyBY.m')

%--------------------------------------------------------------------------
% Chose Fly-by conditions for each pair

Flyby_deltaPumb_and_Cranck=zeros(length(NODEOUT_Matrix(:,1)),2);
Flyby_deltaPumb_and_Cranck(1,:)=[0 0.003];
Flyby_deltaPumb_and_Cranck(2,:)=[0 0.003];
Flyby_deltaPumb_and_Cranck(3,:)=[0 -0.003];
Flyby_deltaPumb_and_Cranck(4,:)=[0 -0.003];
Flyby_deltaPumb_and_Cranck(5,:)=[-0.0010 0.0005];
Flyby_deltaPumb_and_Cranck(6,:)=[+0.0005 +0.0005];
Flyby_deltaPumb_and_Cranck(7,:)=[+0.0005 -0.0005];
Flyby_deltaPumb_and_Cranck(8,:)=[0.000 0.0024];
Flyby_deltaPumb_and_Cranck(9,:)=[0.000 0.0024];
Flyby_deltaPumb_and_Cranck(10,:)=[0.0012 -0.0012];

% Run Refinement of Original Transfer
run script_refine_as_DSMLambert_transfers.m

DV_refined_Original=sum(DVatLeg);

%--------------------------------------------------------------------------
%
NODEIN_Matrix_modified=NODEIN_Matrix;
NODEOUT_Matrix_modified=NODEOUT_Matrix;

for iN=1:length(NODEIN_Matrix_modified(:,1))

    if (iN==1||JumpsTrack(iN)~=0)
        NODEOUT_Matrix_modified(iN,:)=[NODEIN_Matrix_modified(iN,1) ...
            NODEIN_Matrix_modified(iN,2)+Flyby_deltaPumb_and_Cranck(iN,1) ...
            NODEIN_Matrix_modified(iN,3)+Flyby_deltaPumb_and_Cranck(iN,2)];
    else
        NODEOUT_Matrix_modified(iN,:)=[NODEOUT_Matrix_modified(iN-1,1) ...
            NODEOUT_Matrix_modified(iN-1,2)+Flyby_deltaPumb_and_Cranck(iN,1) ...
            NODEOUT_Matrix_modified(iN-1,3)+Flyby_deltaPumb_and_Cranck(iN,2)];
    end


    if (iN==1||JumpsTrack(iN)~=0)
        NODEIN_Matrix_modified(iN,:)=NODEIN_Matrix(iN,:);
    else
        NODEIN_Matrix_modified(iN,:)=NODEOUT_Matrix_modified(iN-1,:);
    end

end


Full_route_new=Full_route;
Full_route_new(indexEnceladus_NODEIN,9:11)=NODEIN_Matrix_modified;
Full_route_new(indexEnceladus_NODEOUT,2:4)=NODEOUT_Matrix_modified(1:end-1,:);

Full_route=Full_route_new;
run script_refine_as_DSMLambert_transfers.m
DV_refined_Modified=sum(DVatLeg);

%% Plotting full GTs

clear Flyby

iNmax=length(NODEIN_Matrix(:,1));
for iN=1:iNmax

NODEIN1=NODEIN_Matrix_modified(iN,:);

%Determine maximum bending due to flyby
rp_flyby  = parsPLOT.INPUTS.Flyby.min_h(1) + parsPLOT.Moon.EquRad(1);           %[km
e_fly     = 1 + ((rp_flyby*NODEIN1(1)^2)/pars.Moon.mu(1));    %[-]
delta_max = 2*asin(1/e_fly);                                      %[rad]
parsPLOT.delta_max = delta_max;


NODEOUT1=NODEOUT_Matrix_modified(iN,:);
disp(['NODE IN : ',num2str(NODEIN1(1)),' km/s ', num2str(NODEIN1(2)),' rads ',num2str(NODEIN1(3)),' rads '])
disp(['NODE OUT : ',num2str(NODEOUT1(1)),' km/s ', num2str(NODEOUT1(2)),' rads ',num2str(NODEOUT1(3)),' rads '])


% Compute flyby parameters
[Flyby(iN)] = Flyby_BuildUp(NODEIN1, NODEOUT1, parsPLOT);
end

colors = cool(length(Flyby));
fig1 = figure( 'Color', [1 1 1] );
hold on;
plotTextureLatLong(parsPLOT.INPUTS.idMoon , parsPLOT.INPUTS.idCentral , 1);
plotSquares(parsPLOT, 1);
axis normal;
for i = 1:size(Flyby, 2)
    Plot_Flyby_GT(Flyby(i), colors(i,:));
end