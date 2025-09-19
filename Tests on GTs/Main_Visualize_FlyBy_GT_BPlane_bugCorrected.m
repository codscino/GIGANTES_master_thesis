%% DEMO SCRIPT TO PLOT GROUNDTRACKS %%
clear all; close all; clc; warning off

% This script should generate and help visualizing a full fly-by with a
% given planet or moon. The idea is to generate both the ground track plot
% and the b-plane trajectory plot.

% obviously the script needs to input some change of node otherwise the
% fly-by is undefined. 

%% DEFINE PLANET, MOON PARAMETERS & CONSTANTS %%
% Define Central Body & Moon of interest
pars.INPUTS.idCentral  = 6; % (5 = Jupiter; 6 = Saturn; 7 = Uranus)
pars.INPUTS.idMoon     = [5]; % (1 = Io; 2 = Europa; 3 = Ganymede; 4 = Callisto)
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
[rrga, vvga] = approxEphem_CC(pars.INPUTS.idMoon, pars.INPUTS.epoch0, pars.INPUTS.idCentral);


% NODEIN Definition
alfain = deg2rad(45); kin = deg2rad(270);
% VinfMust Remain the same and Deflection must be less than pars.delta_max
Delta_Alpha=0.4*delta_max;
Delta_Crank=1.3*delta_max;
Delta_Crank= 0.4356;

% Delta_Alpha=alfain+0.912288518886888;
% Delta_Crank=-kin+5.15955953994486;

NODEIN  = [vinfin, alfain, kin]; %[km/s, rad, rad]

% NODEOut Definition
NODEOUT=[vinfin, alfain+Delta_Alpha, kin+Delta_Crank];

CARNODEIN=vinfAlphaCrank2car(NODEIN, [rrga vvga],  pars.Moon.mu);
CARNODEOUT=vinfAlphaCrank2car(NODEOUT, [rrga vvga],  pars.Moon.mu);
delta = acos( dot(CARNODEIN, CARNODEOUT)./(norm(CARNODEIN)*norm(CARNODEOUT)));


% display info in screen
disp(['NODE IN : ',num2str(NODEIN(1)),' km/s ', num2str(NODEIN(2)),' rads ',num2str(NODEIN(3)),' rads '])
disp(['NODE OUT : ',num2str(NODEOUT(1)),' km/s ', num2str(NODEOUT(2)),' rads ',num2str(NODEOUT(3)),' rads '])
disp(['Deflection angle is ',num2str(delta),' rads'])
disp(['Maximum Deflection angle is ',num2str(delta_max),' rads'])

delta_analytic=acos(cos(NODEIN(2))*cos(NODEOUT(2))+sin(NODEIN(2))*sin(NODEOUT(2))*(cos(Delta_Crank)));

cos_Delta_Crank=(cos(delta_max)-(cos(NODEIN(2))^2*cos(Delta_Alpha)-cos(NODEIN(2))*sin(NODEIN(2))*sin(Delta_Alpha)))/(sin(NODEIN(2))^2*cos(Delta_Alpha)+sin(NODEIN(2))*cos(NODEIN(2))*sin(Delta_Alpha));
Delta_Crank_2max=[acos(cos_Delta_Crank) -acos(cos_Delta_Crank)];

% Compute flyby parameters
[Flyby] = Flyby_BuildUp_BC(NODEIN, NODEOUT, pars);

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
%% PLOT BPLANE DISPLAY OF FLYBY %%

States=Flyby.fly_States;
Distance=sqrt(States(:,1).^2+States(:,2).^2+States(:,3).^2);

MaxDistance=max(Distance); % Define Distance under which the trajectory will be ploted.


% 
% indexDel=find(Distance>MaxDistance);
% States(indexDel,:)=[]; Distance(indexDel)=[]; % deleting states that are above the maximum distance.

[~,indexBpoint]=min(Distance);
Bpoint=States(indexBpoint,1:3);
EntranceHyperbola=States(1,1:3);
ExitHyperbola=States(end,1:3);


%--------------------------------------------------------------------------
% Define The B Plane Components
% The intention is to define a b-plane on which the coordinate Y shows a
% vector within the equator of the moon and perpendicular to the direction
% of motion of the spacecraft
 T=States(1,4:6)/norm(States(1,4:6)); % Tangent Direction (i.e., perpendicular to the B Plane)
 Z=[0 0 1];
 S= crossFast(Z,T)/norm(crossFast(Z,T));
 N =crossFast(T,S);
%--------------------------------------------------------------------------
% Draw Planet: Enceladus
UnitDistance = pars.Moon.OrbRad; % Km, from takubo et al.  
scale = 1;
figure1 = figure('Color',[1 1 1]);
% Create axes
axes1 = axes('Parent',figure1,'ZColor',[1 1 1],'YColor',[1 1 1],...
    'XColor',[1 1 1]);
hold(axes1,'all');
MoonFB='Titan';
[hfig] = drawPlanet(MoonFB,[0 0 0],figure1,scale);
axis equal
view(0,90)
%--------------------------------------------------------------------------
% Draw the BPlane 
% VInfinDir=Flyby.vvinfin/norm(Flyby.vvinfin); %relative direction of motion
CMPos=[0 0 0];
p1=MaxDistance*T;
%    p2=[xEarth-0.3 0 0];
line([ CMPos(1) p1(1)],[ CMPos(2) p1(2)],[ CMPos(3) p1(3)],'Color',[0 0 0])
% Arrow head
% scale factor
H=gca;
xr=get(H,'xlim');
yr=get(H,'ylim');
zr=get(H,'zlim');
s=norm([xr(2)-xr(1),yr(2)-yr(1),zr(2)-zr(1)]);
r=s/300; % radius
h=6*r; % height
n=20; m=20; % grid spacing
[x,y,z]=cylinder(linspace(0,r,n),m);
H=surf(x+p1(1),y+p1(2),z*h+p1(3));
set(H,'FaceColor','k','EdgeColor','k')
P=p1-CMPos; % position vector
U=P/norm(P); % unit vector
D=cross([0,0,1],-U); % rotation axis
if norm(D)<eps, D=[1,0,0]; end
A=acos(dot([0,0,1],-U)); % rotation angle
rotate(H,D,A*180/pi,p1)
%--------------------------------------------------------------------------
CMPos=[0 0 0];
p2=MaxDistance*S;
%    p2=[xEarth-0.3 0 0];
line([ CMPos(1) p2(1)],[ CMPos(2) p2(2)],[ CMPos(3) p2(3)],'Color',[0 0 0])
% Arrow head
% scale factor
H=gca;
xr=get(H,'xlim');
yr=get(H,'ylim');
zr=get(H,'zlim');
s=norm([xr(2)-xr(1),yr(2)-yr(1),zr(2)-zr(1)]);
r=s/300; % radius
h=6*r; % height
n=20; m=20; % grid spacing
[x,y,z]=cylinder(linspace(0,r,n),m);
H=surf(x+p2(1),y+p2(2),z*h+p2(3));
set(H,'FaceColor','k','EdgeColor','k')
P=p2-CMPos; % position vector
U=P/norm(P); % unit vector
D=cross([0,0,1],-U); % rotation axis
if norm(D)<eps, D=[1,0,0]; end
A=acos(dot([0,0,1],-U)); % rotation angle
rotate(H,D,A*180/pi,p2)
%--------------------------------------------------------------------------
CMPos=[0 0 0];
p3=MaxDistance*N;
%    p2=[xEarth-0.3 0 0];
line([ CMPos(1) p3(1)],[ CMPos(2) p3(2)],[ CMPos(3) p3(3)],'Color',[0 0 0])
% Arrow head
% scale factor
H=gca;
xr=get(H,'xlim');
yr=get(H,'ylim');
zr=get(H,'zlim');
s=norm([xr(2)-xr(1),yr(2)-yr(1),zr(2)-zr(1)]);
r=s/300; % radius
h=6*r; % height
n=20; m=20; % grid spacing
[x,y,z]=cylinder(linspace(0,r,n),m);
H=surf(x+p3(1),y+p3(2),z*h+p3(3));
set(H,'FaceColor','k','EdgeColor','k')
P=p3-CMPos; % position vector
U=P/norm(P); % unit vector
D=cross([0,0,1],-U); % rotation axis
if norm(D)<eps, D=[1,0,0]; end
A=acos(dot([0,0,1],-U)); % rotation angle
rotate(H,D,A*180/pi,p3)
%---------------------------------------------------------------------------
% plot fly_by
plot3(States(:,1), States(:,2), States(:,3))
% How periapsis point
line([0 Bpoint(1)],[0 Bpoint(2)], [0 Bpoint(3)])
% Show entrance of the assymptote
line([0 EntranceHyperbola(1)],[0 EntranceHyperbola(2)], [0 EntranceHyperbola(3)])
% Show exit of the assymptote
line([0 ExitHyperbola(1)],[0 ExitHyperbola(2)], [0 ExitHyperbola(3)])
%% Notes

% % Draw the Trajectory plane  first
% XData1=[0 1 1 0];
% YData1=[0 0 1 1];
% ZData1=[0 0 0 0];
% % Create figure
% figure1 = figure;
% % Create axes
% axes1 = axes('Parent',figure1);
% % Create patch
% patch('Parent',axes1,'ZData',ZData1,'YData',YData1,'XData',XData1,...
%     'FaceAlpha',0.1,...
%     'LineWidth',2,...
%     'FaceColor',[1 0 0]);
% view(axes1,[-37.5 30]);

%% Other plots

% Plot distance as a function of time
% figure
% semilogy(Flyby.fly_tts,Distance)

%% Plot the v_infity_sphere of the flyby


INComing_direction=[cos(NODEIN(2)), sin(NODEIN(2))*cos(NODEIN(3)), -sin(NODEIN(2))*sin(NODEIN(3))];
OUTgoing_direction=[cos(NODEOUT(2)), sin(NODEOUT(2))*cos(NODEOUT(3)), -sin(NODEOUT(2))*sin(NODEOUT(3))];

% Create figure
figure1 = figure('WindowState','maximized');

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
quiver3(0,0,0,INComing_direction(1),INComing_direction(2),INComing_direction(3),'LineWidth',2,'Color',[0 0.498039215803146 0],...
    'AutoScaleFactor',0.889999985694885, 'DisplayName','INComing Direction');
quiver3(0,0,0,OUTgoing_direction(1),OUTgoing_direction(2),OUTgoing_direction(3),'LineWidth',2,...
    'Color',[0.850980401039124 0.325490206480026 0.0980392172932625],...
    'AutoScaleFactor',0.889999985694885,'DisplayName','OUTgoing Direction');

PeriapsisPoint=(INComing_direction-OUTgoing_direction)/norm(INComing_direction-OUTgoing_direction);
quiver3(0,0,0,PeriapsisPoint(1),PeriapsisPoint(2),PeriapsisPoint(3),'LineWidth',2,...
     'Color',[0 0 0],...
    'AutoScaleFactor',0.889999985694885,'DisplayName','Periapsis Position');
axis equal
     nfaces = 20;
    [X,Y,Z] = sphere(nfaces);
% 4.2. ploting
surf (X, Y, Z,'Parent',axes1,'FaceColor','none','EdgeColor',[0.800000011920929 0.800000011920929 0.800000011920929]);

RotationDirection=cross(INComing_direction,[0 0 1])/norm(cross(INComing_direction,[0 0 1]));
v1 = eulerAxisAngle(INComing_direction,RotationDirection,delta_max);

thetaRot=linspace(0,2*pi);
SetofAvailableOptions=zeros(length(thetaRot),3);
AvailablePeriapis_alphas_kas=zeros(length(thetaRot),7);
for i=1:length(thetaRot)
    SetofAvailableOptions(i,:)=eulerAxisAngle(v1,INComing_direction,thetaRot(i));
     PossiblePeri=(INComing_direction-SetofAvailableOptions(i,:))/norm(INComing_direction-SetofAvailableOptions(i,:));
     alpha_new1=acos(SetofAvailableOptions(i,1));
     crank_new1=acos(SetofAvailableOptions(i,2)/sin(alpha_new1));
     AvailablePeriapis_alphas_kas(i,:)=[PossiblePeri alpha_new1 crank_new1 -alpha_new1 mod(-crank_new1,2*pi)];
end
plot3(SetofAvailableOptions(:,1),SetofAvailableOptions(:,2),SetofAvailableOptions(:,3),'Parent',axes1,'LineWidth',2,'Color',[0.400000005960464 0 0.400000005960464])

% Create zlabel
zlabel('z (Negative out of plane direction)');
% Create ylabel
ylabel('y (radial outward direction)');
% Create xlabel
xlabel('x (direction of planetary motion)');
view(axes1,[-10.8869329523003 10.4574443141852]);
grid(axes1,'on');
axis(axes1,'tight');
% hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'DataAspectRatio',[1 1 1]);
% Create legend
legend(axes1,'show');



%%
% alfain = deg2rad(45); kin = deg2rad(270);
% Delta_Alpha=0.35*delta_max;
% Delta_Crank=+1.25*delta_max;

% % NODEIN Definition
% alfain = deg2rad(-45); kin = deg2rad(45);
% NODEIN  = [vinfin, alfain, kin]; %[km/s, rad, rad]
% 
% % NODEOut Definition
% % VinfMust Remain the same and Deflection must be less than pars.delta_max
% Delta_Alpha=0.35*delta_max;
% Delta_Crank=1.25*delta_max;

