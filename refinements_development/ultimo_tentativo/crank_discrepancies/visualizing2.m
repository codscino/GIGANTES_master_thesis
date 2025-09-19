% visualize plots with classic GIGANTES script

clc;
clear all;
close all;

pars.GroundTr.t_prop  = 200*60; % minutes, 2000 is more that 1 enceladus full revolution

% pars.EncPlotSize = 1; %true enceladus size, usually too tiny
pars.EncPlotSize = 50;
pars.INPUTS.epoch0  = date2mjd2000([2040 1 8 0 0 0]);

%% ========================================================================
%  1. DEFINE INPUTS & INITIALIZE
%  ========================================================================

% --- Flyby and Central Body Parameters ---
pars.INPUTS.idCentral = 6;      % Central Body: Saturn
pars.INPUTS.idMoon    = 1;      % Enceladus flyby
pars.INPUTS.V_inf     = 4;    % Hyperbolic excess velocity at Enceladus [km/s]
pars.INPUTS.Flyby.min_h = 25;   % Minimum flyby altitude [km]

% % Re-Enter 7:1 Resonance
% nodein  = [4, deg2rad(8.5886), deg2rad(0)];
% nodeout = [4, deg2rad(8.6918), deg2rad(0)];

% nodein =  [4, 0, 0];            % [km/s, rad, rad]
% nodeout = [4, 0.15, deg2rad(1)];   % [km/s, rad, rad]

% Claudio short orbit
% nodein  = [4, deg2rad(45.6918), deg2rad(0)];
% nodeout = [4, deg2rad(45.6918), deg2rad(1)];

% Partial-COT 1 (O/I) - 1st Flyby
nodein  = [4, deg2rad(8.6918), deg2rad(-86.9406)];
nodeout = [4, deg2rad(8.6918), deg2rad(-88.1610)];


% --- Load Gravitational Parameters ---
mu_central_body = getAstroConstants('Saturn', 'Mu');
[~, mu_enceladus, R_enceladus, ~] = satMoonsConstants(1); % Enceladus

% Retrieve Saturn Parameters 
[pars.Planet.mu, pars.Planet.EquRad, pars.Planet.OrbRad, pars.Planet.hmin] = planetConstants(pars.INPUTS.idCentral);

% Retrieve Desired Moon Parameters 
if pars.INPUTS.idCentral == 6
    [pars.Moon.OrbRad, pars.Moon.mu, pars.Moon.EquRad, pars.Moon.hmin] = satMoonsConstants(pars.INPUTS.idMoon);
end
pars.Moon.Vel    = sqrt(pars.Planet.mu/pars.Moon.OrbRad);
pars.Moon.Period = 2*pi*sqrt(pars.Moon.OrbRad^3/pars.Planet.mu);
pars.Moon.HillSph = pars.Moon.OrbRad*( pars.Moon.mu/(3*(pars.Moon.mu + pars.Planet.mu)))^(1/3);

% --- Parameters required for Linked Conic Calculation ---
pars.GroundTr.npoints      = 30000; 

pars.INPUTS.Flyby.hMapping = 300;


%% ========================================================================
%  2. CALCULATE INITIAL PERICENTER STATE FOR PROPAGATION
%  ========================================================================
% [r_enceladus_at_flyby, v_enceladus_at_flyby] = approxEphemSatMoons_cc_numerical(pars.INPUTS.idMoon,  pars.INPUTS.epoch0);
[r_enceladus_at_flyby, v_enceladus_at_flyby] = approxEphem_CC(pars.INPUTS.idMoon,  pars.INPUTS.epoch0, pars.INPUTS.idCentral);

[vvinfin, r0_sc_in, v0_sc_in, ~] = vinfAlphaCrank_to_VinfCART(nodein(1), nodein(2), nodein(3), pars.INPUTS.epoch0, pars.INPUTS.idMoon, pars.INPUTS.idCentral);
[vvinfout, r0_sc_out, v0_sc_out, vvga] = vinfAlphaCrank_to_VinfCART(nodeout(1), nodeout(2), nodeout(3), pars.INPUTS.epoch0, pars.INPUTS.idMoon, pars.INPUTS.idCentral);

rp_flyby  = pars.INPUTS.Flyby.min_h + pars.Moon.EquRad;
e_fly     = 1 + ((rp_flyby*pars.INPUTS.V_inf^2)/pars.Moon.mu);
delta_max = 2*asin(1/e_fly);
pars.delta_max = delta_max;

% b1 = -r0_sc_out./norm(r0_sc_out);
% b3 = cross(r0_sc_out, vvga)./norm(cross(r0_sc_out, vvga));
% b2 = cross(b3,b1);
% Rm = [ b1' b2' b3' ]';

Rm = buildRm(r0_sc_out,vvga); %jose method
% Rm1 = buildRm(r_enceladus_at_flyby, v_enceladus_at_flyby); %my method

[vvinfouBM, delta] = infVelBeforeDefect(vvinfin, vvinfout, delta_max, vvga);
vvinfin_bf = [Rm*vvinfin']';
vvinfouBM_bf = [Rm*vvinfouBM']';

Energy = 0.5*norm(vvinfin_bf)^2;
sma = -mu_enceladus/(2*Energy);
ecc = 1/(sin(delta/2));
rp = sma*(1 - ecc);
hhat = cross(vvinfin_bf, vvinfouBM_bf)./norm(cross(vvinfin_bf, vvinfouBM_bf));
vp = sqrt(norm(vvinfin_bf)^2 + 2*mu_enceladus/rp);

rrp_bf = rp.*(vvinfin_bf - vvinfouBM_bf)./norm(vvinfin_bf - vvinfouBM_bf);
vvp_bf = vp.*cross(hhat, rrp_bf./rp);

Rm_inv = Rm';
rrp_saturn_centric = (Rm_inv * rrp_bf')' + r_enceladus_at_flyby;
vvp_saturn_centric = (Rm_inv * vvp_bf')' + v_enceladus_at_flyby;

rrp_saturn_centric = rrp_saturn_centric(:);
vvp_saturn_centric = vvp_saturn_centric(:);


%% ========================================================================
%  3.1 PROPAGATE THE NBODY TRAJECTORY
%  ========================================================================

duration_sec = pars.GroundTr.t_prop * 60;
time_steps = pars.GroundTr.npoints;
time_vector_fwd = linspace(eps, duration_sec, time_steps)';
time_vector_bwd = linspace(0, -duration_sec, time_steps)';
% 
% ephem_handle = @(t_sec) get_body_positions_wrapper(t_sec, pars.INPUTS.epoch0, actualBodyNaifIDs, spiceParam);
% 
% % Propagate FORWARD for the full time duration
% [time_out_fwd, state_out_fwd, ~, ~, ~] = propagateNBodyODE3(rrp_saturn_centric, vvp_saturn_centric, time_vector_fwd, mu_central_body, mu_TBs, ephem_handle, specialPerturbationIDs, []);
% 
% % Propagate BACKWARD for the full time duration
% [time_out_bwd, state_out_bwd, ~, ~, ~] = propagateNBodyODE3(rrp_saturn_centric, vvp_saturn_centric, time_vector_bwd, mu_central_body, mu_TBs, ephem_handle, specialPerturbationIDs, []);
% 
% % Report propagation completion
% fprintf('Forward propagation completed for %.2f minutes.\n', pars.GroundTr.t_prop);
% fprintf('Backward propagation completed for %.2f minutes.\n', pars.GroundTr.t_prop);
% 
% time_out_bwd = flipud(time_out_bwd);
% state_out_bwd = flipud(state_out_bwd);
% full_time_out_nb = [time_out_bwd; time_out_fwd(1:end)]; %from 2 so the pericentre is not repeated
% full_state_out_nb = [state_out_bwd; state_out_fwd(1:end, :)];
% num_steps = length(full_time_out_nb);

%% ========================================================================
%  3.2 PROPAGATE THE LINKED CONICS TRAJECTORY
%  ========================================================================

% --- Propagate backward from the incoming node for the full time duration ---
[timeLC_bwd, stateLC_bwd, ~, ~, ~] = propagateKeplerODE2(r0_sc_in', v0_sc_in', time_vector_bwd, mu_central_body, []);

% --- Propagate forward from the outgoing node for the full time duration ---
[timeLC_fwd, stateLC_fwd, ~, ~, ~] = propagateKeplerODE2(r0_sc_out', v0_sc_out', time_vector_fwd, mu_central_body, []);

% Flip the backward propagation results
stateLC_bwd_flipped = flipud(stateLC_bwd);
timeLC_bwd_flipped = flipud(timeLC_bwd);

% Merge the two branches using the OUTPUTS from the ODE solver.
full_time_out_LC = [timeLC_bwd_flipped; timeLC_fwd(1:end)];
full_state_out_LC = [stateLC_bwd_flipped; stateLC_fwd(1:end, :)];

num_steps = length(full_state_out_LC);

%% ========================================================================
%  3.3 RESAMPLE TRAJECTORIES TO A COMMON TIME GRID FOR COMPARISON
%  ========================================================================

% The two trajectories have different time steps. To compare them, I must
% resample one onto the time grid of the other. I will use the N-body
% time vector as the "master" reference grid.

% Use interp1 to find the LC state at the times specified by the N-body propagation.
% state_LC_resampled = interp1(full_time_out_LC, full_state_out_LC, full_time_out_nb, 'pchip');

state_LC_resampled = full_state_out_LC;

%% ========================================================================
%  4. PREPARE DATA FOR PLOTTING - INCLUDING TITAN HISTORY
%  ========================================================================

% --- Calculate the full position history of Enceladus ---
r_enc_history = zeros(num_steps, 3);
v_enc_history = zeros(num_steps, 3);
for i = 1:num_steps
    current_mjd = pars.INPUTS.epoch0 + full_state_out_LC(i) / 86400;
    % [r_enc, v_enc] =   approxEphemSatMoons_cc_numerical(pars.INPUTS.idMoon,  current_mjd);
    [r_enc, v_enc] = approxEphem_CC(pars.INPUTS.idMoon,  current_mjd, pars.INPUTS.idCentral);
    r_enc_history(i, :) = r_enc;
    v_enc_history(i, :) = v_enc;
end

% --- Calculate the full position history of Titan ---
r_titan_history = zeros(num_steps, 3);
v_titan_history = zeros(num_steps, 3);
for i = 1:num_steps
    current_mjd = pars.INPUTS.epoch0 + full_state_out_LC(i) / 86400;
    % [r_titan, v_titan] =   approxEphemSatMoons_cc_numerical(5,  current_mjd);
    [r_titan, v_titan] = approxEphem_CC(5,  current_mjd, pars.INPUTS.idCentral);
    r_titan_history(i, :) = r_titan;
    v_titan_history(i, :) = v_titan;
end

fprintf('Ephemeris calculation completed for Enceladus and Titan.\n');


%% ========================================================================
%  9. CREATE ANIMATED PLOTS WITH TITAN AND ENCELADUS
%  ========================================================================

% --- Animated Saturn-Centric Flyby Plot with both Titan and Enceladus ---
liveplot2(state_LC_resampled, r_enc_history, r_titan_history, full_time_out_LC, pars);

% --- Animated Enceladus-Centric Flyby Plot ---
% createAnimatedEnceladusCentricPlot(r_sc_in_flyby_frame_nb, r_sc_in_flyby_frame_LC, r_saturn_in_flyby_frame_nb, full_time_out_nb, pars);



