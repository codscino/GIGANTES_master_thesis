
function [vinfCAR, rr, vv, vvga, kepga] = vinfAlphaCrank_to_VinfCARTClaudio(vinf, alpha, k, t0_mjd, pars)

% % This function computes the v-infinity vector in cartesian coordinates
% % given a set of v-infinity magnitude, pump and crank angles, using spice
% %
% % INPUT
% % - vinf : v-infinity magnitude [km/s]
% % - alpha : pump angle [rad]
% % - k : crank angle [rad]
% % - t0_mjd : epoch of the flyby [MJD2000]
% % - idpl : ID of the flyby body (see also constants.m)
% % - idcentral : ID of the central bodt (see also constants.m)
% %
% % OUTPUT
% % - vinfCAR : 1x3 vector with v-infinity in cartesian coordinates (vinfx,
% % vinfy, vinfz) [km/s]
% % - rr : 1x3 vector of spacecraft position at epoch [km]
% % - vv : 1x3 vector of spacecraft velocity at epoch [km/s]
% % - vvga : 1x3 vector of flyby body velocity at epoch [km/s]
% % - kepga : 1x6 vector of Keplerian elements [a, e, i, omega, Omega, theta]
% %

spiceParam.observer = num2str(pars.INPUTS.NAIFCentral);  % Saturn
% spiceParam.frame = 'IAU_SATURN';
spiceParam.frame    = 'J2000';
spiceParam.abcorr = 'NONE';

% Get moon ephemeris
[rrga, vvga, kepga] = EphSS_car_spice2(pars.INPUTS.NAIFMoon, t0_mjd, true, spiceParam);

% Define v-infinity in TCN frame
% vinfTCN = vinf * [cos(alpha), -sin(alpha)*sin(k), sin(alpha)*cos(k)]; % original
vinfTCN = vinf.*[ sin(alpha)*cos(k), sin(alpha)*sin(k), cos(alpha)];  % Strange pag. 13, eq. 2.28
% vinfTCN = vinf.*[ sin(alpha)*cos(k), cos(alpha), -sin(alpha)*sin(k)];  % Strange pag. 37, eq. 2.156
% vinfTCN = -vinf.*[ sin(alpha)*cos(k), cos(alpha), sin(alpha)*sin(k)]; % Deimos report pag. 25
% vinfTCN = vinf.*[ sin(alpha)*cos(k), -cos(alpha), sin(alpha)*sin(k)]; % Levi, Marazza article pag. 9


% Build proper TCN-to-inertial transformation
% T: Tangential (along velocity)
T_hat = vvga / norm(vvga);

% N: Normal (perpendicular to orbital plane)
C_hat = cross(rrga, vvga);
C_hat = C_hat / norm(C_hat);

% C: Circumferential (completes right-handed system)
N_hat = cross(T_hat, C_hat);

% Transformation matrix from TCN to inertial
% TCN_to_inertial = [T_hat', C_hat', N_hat'];
TCN_to_inertial = [N_hat', C_hat', T_hat'];

% Transform v-infinity to inertial frame
vinfCAR = (TCN_to_inertial * vinfTCN')';

% Output spacecraft state
rr = rrga;
vv = vvga + vinfCAR;

end