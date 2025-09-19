% script_refine_as_DSMLambert_transfers.m

nLegs=length(Full_route(:,1));

DVatLeg=zeros(1,nLegs);
for iLeg=1:nLegs

TOF_seconds = Full_route(iLeg,7)*3600*24; %[seconds]

%Departure Position & Velocity
DepFlag=Full_route(iLeg,1);


if DepFlag==5
%    Titan position & velocity
    MAnomaly_MoonDeparture=Full_route(iLeg,5);
    SVatTitan=SVTitan_func(MAnomaly_MoonDeparture);
    r_moon_departure=SVatTitan(1:3);
    v_moon_departure=SVatTitan(4:6);
elseif DepFlag==1
%    Enceladus position & velocity
    MAnomaly_MoonDeparture=Full_route(iLeg,6);
    SVatEnceladus=SVEnceladus_func(MAnomaly_MoonDeparture);
    r_moon_departure=SVatEnceladus(1:3);
    v_moon_departure=SVatEnceladus(4:6);
end


%Arrival Position & Velocity
ArrFlag=Full_route(iLeg,8);
if ArrFlag==5
%    Titan position & velocity
    MAnomaly_MoonArrival=Full_route(iLeg,12);
    SVatTitan=SVTitan_func(MAnomaly_MoonArrival);
    r_moon_arrival=SVatTitan(1:3);
    v_moon_arrival=SVatTitan(4:6);
elseif ArrFlag==1
%    Enceladus position & velocity
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


% [dv]=Wrapper_GALambertDSM_refinement(MOONDepArrIndex,NODEOUT1, sv1_planet, NODEOUT2, sv2_vector, TOF_seconds, 0.01, pars)

disp(['Optimization of leg ',num2str(iLeg),' out of ',num2str(nLegs)])
%--------------------------------------------------------------------------
%Step 1 - Define Limits of Design Variables
LB =0; UB =1;
%--------------------------------------------------------------------------
%Step 1- Create quick full exploration.
fun=@(x)Wrapper_GALambertDSM_refinement(MOONDepArrIndex,NODEOUT1, sv1_planet, NODEOUT2, sv2_vector, TOF_seconds, x(1), pars);
%--------------------------------------------------------------------------
%step 3- optimize
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
% 


DVatLeg(iLeg)=Out_DSM_GA;


end

sum(DVatLeg)