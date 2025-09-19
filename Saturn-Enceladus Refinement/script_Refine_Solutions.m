 clear all; close all; clc; format short g; warning off;
% % script_Refine_Solutions.m
% %--------------------------------------------------------------------------
% load('Saturn_Enceladus_Titan_exploration_results_full3.mat')
% %%
% numEnceladusFB=zeros(1,length(treeExploration));
% numTitanFB=zeros(1,length(treeExploration));
% Delta_V=zeros(1,length(treeExploration));
% for iN=1:length(treeExploration)
%     % Analyse each individual route
%     Entire_route=treeExploration{iN}.route;
% 
%     % Enceladus FlyBys
%     indexRelevantFlyBys=find(Entire_route(:,8)==1);
%     numEnceladusFB(iN)=length(indexRelevantFlyBys);
% 
%     % Overall DV Cost of the Route
%     routeDVCosts=treeExploration{iN}.routeDVCosts;
%     routeDVCosts(find(isnan(routeDVCosts)))=[];
%     Delta_V(iN)=sum(routeDVCosts);
% 
%     % Titan FlyBys
%     indexRelevantFlyBys=find(Entire_route(:,8)==5);
%     numTitanFB(iN)=length(indexRelevantFlyBys);
% %--------------------------------------------------------------------------
% % Enhance TreeExploration Structure
% treeExploration{iN}.numEnceladusFB=numEnceladusFB(iN);
% treeExploration{iN}.TotalDVCost=Delta_V(iN);
% end
% nFBEnceladusMin=6; MaxDV_Allowed=0.1;nFBTitanMin=2;
% % remove cases with only one Titan encounter (departure)
% Good4Titan=find(numTitanFB>=nFBTitanMin);
% treeExploration = treeExploration(Good4Titan);
% numEnceladusFB=numEnceladusFB(Good4Titan);
% Delta_V=Delta_V(Good4Titan);
% % remove cases with low number of enceladus flybys
% Good4Enceladus=find(numEnceladusFB>=nFBEnceladusMin);
% treeExploration = treeExploration(Good4Enceladus);
% numEnceladusFB=numEnceladusFB(Good4Enceladus);
% Delta_V=Delta_V(Good4Enceladus);
% % Remove Cases with high DB
% undexMaxDV=find(Delta_V<MaxDV_Allowed);
% treeExploration = treeExploration(undexMaxDV);
% numEnceladusFB=numEnceladusFB(undexMaxDV);
% Delta_V=Delta_V(undexMaxDV);
% %--------------------------------------------------------------------------
% %% Postprocessing and ranking
% Score_ScientificInteres=zeros(1,length(treeExploration));
% Delta_longitude=zeros(1,length(treeExploration));
% Total_time_of_flight=zeros(1,length(treeExploration));
% MaxInclination_visited=zeros(1,length(treeExploration));
% for iN=1:length(treeExploration)
%     % Analyse each individual route
%     Entire_route=treeExploration{iN}.route;
% 
%     % Enceladus FlyBys
%     indexRelevantFlyBys=find(Entire_route(:,8)==1);
% 
%     % Diversity on Flybys - Currently diversity is considered only by
%     % different alpha and cranks. The idea is that different speeds of approach  won't
%     % change too much the fly-by, instead different angles of approach
%     % will.
%     RelevantNodes=[Entire_route(indexRelevantFlyBys,9) Entire_route(indexRelevantFlyBys,10) Entire_route(indexRelevantFlyBys,11)];
%     Unique_nodes = table(round(RelevantNodes(:,2)/RAD/5),round(RelevantNodes(:,3)/RAD/5)); 
%     Score_ScientificInteres(iN)=height(unique(Unique_nodes));
% 
%     % Longitude of furthest Enceladus FBs - relevant to illumination
%     % Conditions
%     Initial_Longitude=Entire_route(indexRelevantFlyBys(1),13);
%     LongitudeofallEnceladusNodes=Entire_route(indexRelevantFlyBys,13);
%     Delta_longitude(iN)=max(abs((Initial_Longitude-LongitudeofallEnceladusNodes)))/RAD;
% 
% 
%     %Overall Time of Flight
%     Total_time_of_flight(iN)=sum(Entire_route(:,7)); % days
% 
%     % maximum inclination attained
%     vEnceladus=sqrt(pars.Planet.mu/pars.Moon.OrbRad(1));
%     SVMoon_naive=[pars.Moon.OrbRad(1) 0 0 0 vEnceladus 0];
%     inclination_deg=zeros(1,length(RelevantNodes(:,1)));
%     for jN=1:length(RelevantNodes(:,1))
%     [vinfCAR_TITAN_OUT, rr, vv] = vinfAlphaCrank2car(RelevantNodes(1,:), SVMoon_naive, pars.Planet.mu); 
%     ElementsOrbitSc=car2kep([rr, vv],pars.Planet.mu);
%     inclination_deg(jN)=ElementsOrbitSc(3)/RAD;
%     end
%     MaxInclination_visited(iN)=max(inclination_deg);
% %--------------------------------------------------------------------------
% % Enhance TreeExploration Structure
% treeExploration{iN}.Score_ScientificInteres= Score_ScientificInteres(iN);
% treeExploration{iN}.Delta_longitude_excursion=Delta_longitude(iN);
% treeExploration{iN}.TotalTimeofFlight=Total_time_of_flight(iN);
% end
% %--------------------------------------------------------------------------
% %%
% indexMax_Score=find(Score_ScientificInteres==max(Score_ScientificInteres));
% 
% treeExploration_SubSet=treeExploration(indexMax_Score);
% indexValue=[1:1:length(indexMax_Score)];
% InfoSubSet=[indexValue' Delta_V(indexMax_Score)' Total_time_of_flight(indexMax_Score)' Delta_longitude(indexMax_Score)']
% 
% 
% [~,indexSorted]=sort(InfoSubSet(:,3));
% 
% InfoSubSet=InfoSubSet(indexSorted,:);
% 
% 
% % ploting the maximum Diversity  tour
% IndexCheck=1705
% AFull_route=treeExploration_SubSet{IndexCheck}.route;
% totalDV=sum(treeExploration_SubSet{IndexCheck}.routeDVCosts)
% routeAnatomy=treeExploration_SubSet{IndexCheck}.routeAnatomy;
% plotfullSaturnTransfer(AFull_route,pars)

load('all.mat')

%% Refinement 

% this fuction must rereun the precomputed full route. Note the end node
% and starting node on the same planet may now have a small discontinutity
% due the fly-by itself. 
% since Lambert arc won't work, a MGA+DSM model is used, where one variable
% is optimized. This will be the fraction of time to the DSM, since
% everything else is already prefixed. 

Transfer2Refine.Full_route=AFull_route;

% First step it is to reconctruct the DV and full transfer. 
%--------------------------------------------------------------------------
% Begining of the function here
%--------------------------------------------------------------------------
% Constants & Parameters
RAD=pi/180;
Days2Sec=3600*24;
nEnceladus=sqrt(pars.Planet.mu/pars.Moon.OrbRad(1)^3);
nTitan=sqrt(pars.Planet.mu/pars.Moon.OrbRad(2)^3);
%--------------------------------------------------------------------------
% Epremeris Moons Set up
KepTitan=[pars.Moon.OrbRad(2) 0 0 0 0 0];
KepEnceladus=[pars.Moon.OrbRad(1) 0 0 0 0 0];
SVTitan_func = @(TrA)kep2car([KepTitan(1:5) TrA],pars.Planet.mu);
SVEnceladus_func = @(TrA)kep2car([KepEnceladus(1:5) TrA],pars.Planet.mu);
%--------------------------------------------------------------------------


Full_route=Transfer2Refine.Full_route;


% indexEnceladus=find(Full_route(:,1)==1);
% 
% 
% Full_route(2,3)=Full_route(2,3)+0.0015; 
% Full_route(5,4)=Full_route(5,4)+0.0015;
% Full_route(5,9:11)=Full_route(5,2:4);
% Full_route(6,4)=Full_route(6,4)+0.003;
% Full_route(6,9:11)=Full_route(6,2:4);
% Full_route(7,4)=Full_route(7,4)+0.0045;
% Full_route(7,9:11)=Full_route(7,2:4);
% Full_route(8,4)=Full_route(8,4)+0.003;
% Full_route(8,9:11)=Full_route(8,2:4);
% Full_route(9,4)=Full_route(9,4)+0.0015;
% Full_route(9,9:11)=Full_route(9,2:4);
% Full_route(10,4)=Full_route(10,4)+0.000;



nLegs=length(Full_route(:,1));

DVatLeg=zeros(1,nLegs);
for iLeg=1:nLegs

TOF_seconds = Full_route(iLeg,7)*3600*24; %[seconds]

% Departure Position & Velocity
DepFlag=Full_route(iLeg,1);


if DepFlag==5
    % Titan position & velocity
    MAnomaly_MoonDeparture=Full_route(iLeg,5);
    SVatTitan=SVTitan_func(MAnomaly_MoonDeparture);
    r_moon_departure=SVatTitan(1:3);
    v_moon_departure=SVatTitan(4:6);
elseif DepFlag==1
    % Enceladus position & velocity
    MAnomaly_MoonDeparture=Full_route(iLeg,6);
    SVatEnceladus=SVEnceladus_func(MAnomaly_MoonDeparture);
    r_moon_departure=SVatEnceladus(1:3);
    v_moon_departure=SVatEnceladus(4:6);
end


% Arrival Position & Velocity
ArrFlag=Full_route(iLeg,8);
if ArrFlag==5
    % Titan position & velocity
    MAnomaly_MoonArrival=Full_route(iLeg,12);
    SVatTitan=SVTitan_func(MAnomaly_MoonArrival);
    r_moon_arrival=SVatTitan(1:3);
    v_moon_arrival=SVatTitan(4:6);
elseif ArrFlag==1
    % Enceladus position & velocity
    MAnomaly_MoonArrival=Full_route(iLeg,13);
    SVatEnceladus=SVEnceladus_func(MAnomaly_MoonArrival);
    r_moon_arrival=SVatEnceladus(1:3);
    v_moon_arrival=SVatEnceladus(4:6);
end



% Alternative Calculation

NODEOUT1=Full_route(iLeg,2:4);
sv1_planet=[r_moon_departure v_moon_departure];

if iLeg==nLegs
NODEOUT2=Full_route(iLeg,9:11);   
else
NODEOUT2=Full_route(iLeg+1,2:4);
end

sv2_vector=[r_moon_arrival v_moon_arrival];
MOONDepArrIndex=[DepFlag ArrFlag];



[dv]=Wrapper_GALambertDSM_refinement(MOONDepArrIndex,NODEOUT1, sv1_planet, NODEOUT2, sv2_vector, TOF_seconds, 0.01, pars)

disp(['Optimization of leg ',num2str(iLeg),' out of ',num2str(nLegs)])
%--------------------------------------------------------------------------
% Step 1 - Define Limits of Design Variables
LB =0; UB =1;
%--------------------------------------------------------------------------
% Step 1- Create quick full exploration.
fun=@(x)Wrapper_GALambertDSM_refinement(MOONDepArrIndex,NODEOUT1, sv1_planet, NODEOUT2, sv2_vector, TOF_seconds, x(1), pars);
%--------------------------------------------------------------------------
% step 3- optimize
GenLimit = 500;
StallGenLimit_VAL = 20;
optionsga = optimoptions(@ga,...
    'Generations', GenLimit,...
    'PopulationSize',500,...
    'TolFun', 1e-4,...
    'Display','off',...
    'StallGenLimit',StallGenLimit_VAL);
nres=length(LB);
[xSol_DSM,Out_DSM_GA,exitflag] = ga(@(x)fun(x),nres,[],[],[],[],LB,UB,[],[],optionsga);
disp(['Leg ',num2str(iLeg),' out of ',num2str(nLegs),' requires a DV of ',num2str(Out_DSM_GA),' km/s'])
%% 


DVatLeg(iLeg)=Out_DSM_GA;


end

sum(DVatLeg)




%% Plotting full GTs
% Define Central Body & Moon of interest
parsPLOT.INPUTS.idCentral  = 6; % (5 = Jupiter; 6 = Saturn; 7 = Uranus)
parsPLOT.INPUTS.idMoon     = 1; % (1 = Io; 2 = Europa; 3 = Ganymede; 4 = Callisto)
% (1 = Enceladus; 2 = Thetys; 3 = Dione; 4 = Rhea; 5 = Titan)
% (1 = Miranda; 2 = Ariel; 3 = Umbriel; 4 = Titania; 5 = Oberon)

[parsPLOT.Planet.mu, parsPLOT.Planet.EquRad, parsPLOT.Planet.OrbRad, parsPLOT.Planet.hmin] = planetConstants(parsPLOT.INPUTS.idCentral); %[km3/s2],[km],[km] & [km]

% Retrieve Desired Moon Parameters
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

%--------------------------------------------------------------------------
% Load some necessary parameters
% Define parameters regarding the flyby
parsPLOT.INPUTS.Flyby.min_h    = 25;   %[km] Minimum flyby altitude
parsPLOT.GroundTr.npoints      = 30e3; % N° of time steps for ground-track propagation
parsPLOT.GroundTr.t_prop       = 5;    %[minutes] Time of flyby hyperbola propagation
parsPLOT.INPUTS.Flyby.hMapping = 300;  %[km] Max altitude to consider mapping

% Define starting epoch
parsPLOT.INPUTS.epoch0 = 0;

% Load sectirs for groundtrack plots
load('sectorObj.mat');
load('sectorObj_.mat');
parsPLOT.sectorObj  = sectorObj;
parsPLOT.sectorObj_ = sectorObj_;
parsPLOT.rect       = defineRectangleMapping(); % --> save the rectangle mapping

%--------------------------------------------------------------------------
% Step 1. Identify all the In and out nodes at Enceladus
% 1.1 In nodes at Enceladus
indexINNODES=find(Full_route(:,8)==1);
INNODES_LIST=Full_route(indexINNODES,9:11);
% 1.2. Out nodes at Enceladus
indexOUTNODES=find(Full_route(:,1)==1);
OUTNODES_LIST=Full_route(indexOUTNODES,2:4);
%--------------------------------------------------------------------------
iNmax=length(indexOUTNODES);
for iN=1:iNmax

NODEIN1=INNODES_LIST(iN,:);

% Determine maximum bending due to flyby
rp_flyby  = parsPLOT.INPUTS.Flyby.min_h(1) + parsPLOT.Moon.EquRad(1);           %[km
e_fly     = 1 + ((rp_flyby*NODEIN1(1)^2)/pars.Moon.mu(1));    %[-]
delta_max = 2*asin(1/e_fly);                                      %[rad]
parsPLOT.delta_max = delta_max;


NODEOUT1=OUTNODES_LIST(iN,:);
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