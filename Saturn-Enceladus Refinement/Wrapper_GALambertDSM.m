function [DV_Total]=Wrapper_GALambertDSM(MOONDepArrIndex,NODEIN1,NODEOUT1, sv1_planet, NODEOUT2, sv2_vector, ToF, alpha, pars)

r2_planet=sv2_vector(1:3);

% first issue in the function is that it does not consider multi
% revolution.
[Dsm_vector,v2]=GALambertDSM(NODEOUT1,sv1_planet, r2_planet, ToF, alpha, pars.Planet.mu);
%--------------------------------------------------------------------------
% DV Defect at Enceladus

[vinfCAR_In] = vinfAlphaCrank2car(NODEIN1, sv1_planet, pars.Planet.mu);

vinf_minus_departure=vinfCAR_In;

[vinfCAR] = vinfAlphaCrank2car(NODEOUT1, sv1_planet, pars.Planet.mu);
vinf_plus_departure=vinfCAR;

DepMoon=MOONDepArrIndex(1);
if DepMoon==5; DepMoon=2; end
rp_min=pars.Moon.EquRad(DepMoon)+pars.Moon.hmin(DepMoon);
ecchyperbolic=1+rp_min*norm(vinf_minus_departure)^2/pars.Moon.mu(DepMoon);
DeltaMax_Dep=2*asin(1/ecchyperbolic); % maximum deflection achieved during a flyby
% necessary deflection
delta_at_Dep=acos(dot(vinf_minus_departure,vinf_plus_departure)/norm(vinf_minus_departure)/norm(vinf_plus_departure));
if delta_at_Dep<=DeltaMax_Dep
    DV_at_Dep=abs(norm(vinf_plus_departure)-norm(vinf_minus_departure));
else
    DV_at_Dep=sqrt(norm(vinf_plus_departure)^2+norm(vinf_minus_departure)^2-2*norm(vinf_plus_departure)*norm(vinf_minus_departure)*cos(DeltaMax_Dep-delta_at_Dep));
end
%--------------------------------------------------------------------------
% DV DSM
DV_at_DSM=norm(Dsm_vector);
%--------------------------------------------------------------------------
% DV at Titan
% Calculation of dv_defect for closest match possible (NOTE it may not be necessarily feasible)

v_inf_plus_titan  = vinfAlphaCrank2car(NODEOUT2, sv2_vector, pars.Planet.mu);
v_inf_minus_titan = v2-sv2_vector(4:6);
%--------------------------------------------------------------------------
% Maximum Deflection for final Fly-by
ArrMoon=MOONDepArrIndex(2);
if ArrMoon==5; ArrMoon=2; end
rp_min=pars.Moon.EquRad(ArrMoon)+pars.Moon.hmin(ArrMoon);
ecchyperbolic=1+rp_min*norm(v_inf_minus_titan)^2/pars.Moon.mu(ArrMoon);
DeltaMax_Titan=2*asin(1/ecchyperbolic); % maximum deflection achieved during a flyby
% necessary deflection
delta_at_Titan=acos(dot(v_inf_minus_titan,v_inf_plus_titan)/norm(v_inf_minus_titan)/norm(v_inf_plus_titan));
if delta_at_Titan<=DeltaMax_Titan
    DV_at_Arr=abs(norm(v_inf_plus_titan)-norm(v_inf_minus_titan));
else
    DV_at_Arr=sqrt(norm(v_inf_plus_titan)^2+norm(v_inf_minus_titan)^2-2*norm(v_inf_plus_titan)*norm(v_inf_minus_titan)*cos(DeltaMax_Titan-delta_at_Titan));
end
%--------------------------------------------------------------------------
% Total DV of Transfer
DV_Total=DV_at_Dep+DV_at_DSM+DV_at_Arr;
end
%--------------------------------------------------------------------------
function [Dsm_vector,v2,LambertCases]=GALambertDSM(NODE,sv1_planet, r2_planet, ToF, alpha, muPlanet)

% 1. Compute State Vector at departure
[~, rr, vv] = vinfAlphaCrank2car(NODE, sv1_planet, muPlanet);
sv_out=[rr vv];
% 2. Time of flight of propagation
ToF2DSM=alpha*ToF;
% 3. Propagate to Lambert Initial Point.
[outDSM,errorflag] = kepPro(sv_out,ToF2DSM, muPlanet);
if errorflag==1
    error
end

LambertCases=[0 0; 1 0; 1 1; 2 0; 2 1];
nLC=length(LambertCases(:,1));
v1=zeros(nLC,3); v2=zeros(nLC,3);
Dsm_vector=zeros(nLC,3);
Dsm_norm=zeros(nLC,1);
for iLC=1:nLC
    % 4.  Compute Lambert arc
    orbitType=0;
    Nrev=LambertCases(iLC,1);
    Ncase=LambertCases(iLC,2);
    optionsLMR=0;
    [~,~,~,ERROR,v1(iLC,:),v2(iLC,:)] = lambertMR(outDSM(1:3),r2_planet,(1-alpha)*ToF,muPlanet,orbitType,Nrev,Ncase,optionsLMR);
    %--------------------------------------------------------------------------
    Dsm_vector(iLC,:)=v1(iLC,:)-outDSM(4:6);
    Dsm_norm(iLC)=norm(Dsm_vector(iLC,:));
end
[minDV,indexMin]=min(Dsm_norm); 

Dsm_vector=Dsm_vector(indexMin,:);
v2=v2(indexMin,:);
LambertCases=LambertCases(indexMin,:);

if ERROR~=0
    warning('LAMBERT DID NOT CONVERGE')
end

end