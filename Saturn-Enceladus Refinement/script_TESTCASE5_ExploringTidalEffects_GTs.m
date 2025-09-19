clear all; close all; clc; format short g; warning off;
%% TESTCASE 5
%--------------------------------------------------------------------------
% Here I search for 10 Ecneladus flybys, all equatorial or close to
% equatorial and at different true longitudes simply to explore different
% tidal configurations. 

% Hence I search for 
%%
%--------------------------------------------------------------------------
load('Saturn_Enceladus_Titan_exploration_results_full5.mat')
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
%%
numEnceladusFB=zeros(1,length(treeExploration));
Delta_V=zeros(1,length(treeExploration));
for iN=1:length(treeExploration)
    % Analyse each individual route
    Entire_route=treeExploration{iN}.route;

    % Enceladus FlyBys
    indexRelevantFlyBys=find(Entire_route(:,8)==1);
    numEnceladusFB(iN)=length(indexRelevantFlyBys);

    % Overall DV Cost of the Route
    Delta_V(iN)=treeExploration{iN}.TotalDVCost;
    %--------------------------------------------------------------------------
    % Enhance TreeExploration Structure
    treeExploration{iN}.numEnceladusFB=numEnceladusFB(iN);
end
nFBEnceladusMin=3; MaxDV_Allowed=0.1;
% remove cases with low number of enceladus flybys
Good4Enceladus=find(numEnceladusFB>=nFBEnceladusMin);
treeExploration = treeExploration(Good4Enceladus);
numEnceladusFB=numEnceladusFB(Good4Enceladus);
Delta_V=Delta_V(Good4Enceladus);
%% New filtering
Score_ScientificInteres=zeros(1,length(treeExploration));
Delta_longitude=zeros(1,length(treeExploration));
Total_time_of_flight=zeros(1,length(treeExploration));
MaxInclination_visited=zeros(1,length(treeExploration));
Score_TidalSpread=zeros(1,length(treeExploration));
for iN=1:length(treeExploration)
    % Analyse each individual route
    Entire_route=treeExploration{iN}.route;

    % Enceladus FlyBys
    indexRelevantFlyBys=find(Entire_route(:,8)==1);

    % Diversity on Flybys - Currently diversity is considered only by
    % different alpha and cranks. The idea is that different speeds of approach  won't
    % change too much the fly-by, instead different angles of approach
    % will.
    RelevantNodes=[Entire_route(indexRelevantFlyBys,9) Entire_route(indexRelevantFlyBys,10) Entire_route(indexRelevantFlyBys,11)];
    Unique_nodes = table(round(RelevantNodes(:,2)/RAD/5),round(RelevantNodes(:,3)/RAD/5));
    Score_ScientificInteres(iN)=height(unique(Unique_nodes));

    % Longitude of furthest Enceladus FBs - relevant to illumination
    % Conditions
    Initial_Longitude=Entire_route(indexRelevantFlyBys(1),13);
    LongitudeofallEnceladusNodes=Entire_route(indexRelevantFlyBys,13);
    Delta_longitude(iN)=max(abs((Initial_Longitude-LongitudeofallEnceladusNodes)))/RAD;


    %Overall Time of Flight
    Total_time_of_flight(iN)=sum(Entire_route(:,7)); % days

    % True Longitude of FlyBys
    LongitudeofallEnceladusNodes=Entire_route(indexRelevantFlyBys,13);
    Unique_nodes = table(round(LongitudeofallEnceladusNodes/RAD/15));
    Score_TidalSpread(iN)=height(unique(Unique_nodes));

    Travel_matrix=zeros(length(LongitudeofallEnceladusNodes)); 
    for itm=1:length(LongitudeofallEnceladusNodes)
        for jtm=itm:length(LongitudeofallEnceladusNodes)
            Travel_matrix(itm,jtm)=abs(LongitudeofallEnceladusNodes(itm)-LongitudeofallEnceladusNodes(jtm));
        end
    end
    Score_TidalSpread2(iN)=height(unique(Unique_nodes))*sum(Travel_matrix(:));

    % maximum inclination attained
    vEnceladus=sqrt(pars.Planet.mu/pars.Moon.OrbRad(1));
    SVMoon_naive=[pars.Moon.OrbRad(1) 0 0 0 vEnceladus 0];
    inclination_deg=zeros(1,length(RelevantNodes(:,1)));
    for jN=1:length(RelevantNodes(:,1))
        [vinfCAR_TITAN_OUT, rr, vv] = vinfAlphaCrank2car(RelevantNodes(1,:), SVMoon_naive, pars.Planet.mu);
        ElementsOrbitSc=car2kep([rr, vv],pars.Planet.mu);
        inclination_deg(jN)=ElementsOrbitSc(3)/RAD;
    end
    MaxInclination_visited(iN)=max(inclination_deg);
    %--------------------------------------------------------------------------
    % Enhance TreeExploration Structure
    treeExploration{iN}.Score_ScientificInteres= Score_ScientificInteres(iN);
    treeExploration{iN}.Delta_longitude_excursion=Delta_longitude(iN);
    treeExploration{iN}.TotalTimeofFlight=Total_time_of_flight(iN);
end
%--------------------------------------------------------------------------
%%

% Create figure
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
% Create plot
plot(Score_TidalSpread,numEnceladusFB,'MarkerFaceColor',[0 0 0],'Marker','o','LineStyle','none');
% Create ylabel
ylabel('\Deltav [km/s]');
% Create xlabel
xlabel('\DeltaL [degrees]');
box(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',12,'FontWeight','bold',...
    'XGrid','on','YGrid','on');

% ploting the maximum longitude tour
[MaxL,indexMax]=max(Score_TidalSpread);
Full_route=treeExploration{indexMax}.route;
totalDV=sum(treeExploration{indexMax}.routeDVCosts)
routeAnatomy=treeExploration{indexMax}.routeAnatomy;
plotfullSaturnTransfer(Full_route,pars)




indexMax_Score=find(Score_TidalSpread==4); % we seek for low DV transfers with only one type of fly-by. 

treeExploration_SubSet=treeExploration(indexMax_Score);
indexValue=[1:1:length(indexMax_Score)];
InfoSubSet=[indexValue' numEnceladusFB(indexMax_Score)' Delta_V(indexMax_Score)'  Total_time_of_flight(indexMax_Score)' Score_TidalSpread2(indexMax_Score)']


[~,indexSorted]=sort(InfoSubSet(:,5),'descend');

InfoSubSet=InfoSubSet(indexSorted,:);

% ploting the maximum Diversity  tour
IndexCheck=29277 % This is the exemple to reanalize.
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

 indexNode=10;
 
run('Plot_singleFlyBY.m')

NODEIN  = NODEIN_Matrix(indexNode,:); %[km/s, rad, rad]
% display info in screen
% Determine maximum bending due to flyby
rp_flyby  = parsPLOT.INPUTS.Flyby.min_h + parsPLOT.Moon.EquRad;           %[km
e_fly     = 1 + ((rp_flyby*NODEIN(1)^2)/parsPLOT.Moon.mu);    %[-]
delta_max = 2*asin(1/e_fly);                                      %[rad]
parsPLOT.delta_max = delta_max;
%--------------------------------------------------------------------------
% Chose Fly-by conditions for each pair

Flyby_deltaPumb_and_Cranck=zeros(length(NODEOUT_Matrix(:,1)),2);
Flyby_deltaPumb_and_Cranck(1,:)=[0.00297537806887696 -0.00462388698012824];
Flyby_deltaPumb_and_Cranck(2,:)=[-0.00297537806887696 -0.00462388698012824];
Flyby_deltaPumb_and_Cranck(3,:)=[0.00316733794428838 6.66400187462506e-08];
Flyby_deltaPumb_and_Cranck(4,:)=[-0.00104575531788651 -0.000490465429870121];
Flyby_deltaPumb_and_Cranck(5,:)=[-0.00104575531788651 0.000490465429870121];
Flyby_deltaPumb_and_Cranck(6,:)=[0.00181631186790815 0];
Flyby_deltaPumb_and_Cranck(7,:)=[-0.000799659341213439 -0.00150501325980350];
Flyby_deltaPumb_and_Cranck(8,:)=[0.000799659341213439 0.00150501325980350];
Flyby_deltaPumb_and_Cranck(9,:)=[0.00146352747354158 0.000356093008840613];
Flyby_deltaPumb_and_Cranck(10,:)=[0.000919111256327905 0.00131177066646001];
Flyby_deltaPumb_and_Cranck(11,:)=[-0.000919111256327905 -0.00131252096238307];
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