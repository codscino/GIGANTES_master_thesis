function F = roots_Vinf_Flipping(x, Reso, N_fly, pars)

% Equations 1 and 2 are a system of nonlinear equations which have to be
% solved to retrieve the Vinf & pump angle associated to the orbit which
% allows to perform a V-infinity "flip" (i.e: change the spacecraft orbit
% from Inbound to Outbound or viceversa with a single flyby). The equations
% being solved are based on Buffington et al. - Global Moon Coverage via
% Hyperbolic Flybys, 23rd International Symposium on Space Flight Dynamics.

% Author: Jose Carlos Garcia Mateas
% Last revision: 24/06/2024

%% INPUTS %%
% x: Variables of the nonlinear equations to be solved
%    x(1) = Vinf [km/s], x(2) = Pump angle [rad]
%
% Reso: Resonance ratio, where Reso(1) = N° of Moon orbital periods &
%       Reso(2) = N° of spacecraft orbital periods
%
% N_fly: N° of flybys to be performed by the spacecraft
%
% pars: Structure containing problem parameters & constants

%% CHANGELOG %%
% 24/06/2024, J.C Garcia Mateas: revised function & added comments.

%% FUNCTION %%

F(1) = (1 + (pars.rp_flyby*x(1)^2)/pars.Moon.mu)*sin(x(2)) - sin(pi/(2*N_fly));

F(2) = (pars.Moon.Vel^2*(2 - (Reso(2)/Reso(1))^(2/3)) - x(1)^2 - pars.Moon.Vel^2)/(2*x(1)*pars.Moon.Vel) - cos(x(2));

end