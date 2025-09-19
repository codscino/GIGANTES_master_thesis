function [RA, DEC, W] = enceladus_orientation(t)
% enceladus_orientation Calculates Enceladus orientation angles (RA, DEC, W)
%
% This function computes the inertial orientation of Enceladus (spin-pole
% right ascension RA spin-pole DEC, and prime meridian
% angle W) in the International Celestial Reference Frame (ICRF) based on
% the model described in Park et al. (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023JE008054)
% the model is also described in a short form on a spice PCK kernel 
% https://naif.jpl.nasa.gov/pub/naif/CASSINI/kernels/pck/enceladus_ssd_230702_v1.tpc
%
% INPUT:
%   t : days in MJD2000
%
% OUTPUT:
%   RA     : Right Ascension of the spin pole (degrees)
%   DEC    : Declination of the spin pole (degrees)
%   W      : Prime Meridian angle (rotation around spin pole) (degrees)
%
% The formulas from the paper:
% RA(t) = RA0 + dRA*T + sum(Aj*sin(Sj)) for j=1 to 17
% DEC(t) = DEC0 + dDEC*T + sum(Bj*cos(Sj)) for j=1 to 17
% W(t) = W0  + dW*d + sum(Cj*sin(Sj)) for j=1 to 18
% where:
%   d = time in days since J2000.0
%   T = time in Julian centuries since J2000.0 (1 Julian century = 36525 days)
%   Sj are nutation arguments of the form Sj0 + Sj_rate * T

% --- Time conversion ---
d = t; % Days since J2000.0
T = d / 36525; % Julian centuries since J2000.0

% --- Constants for linear terms (from paper equations / SPICE POLE values) ---
RA0 = 40.592915;
dRA = -0.0902111773;

DEC0 = 83.534180;
dDEC = -0.0071054901;

 % W0 is measured from -5E longitude crater Salih(Park et al.)
W0 = 7.120600; % this values is from spice, it can change in the future with more observations
% W0 = W0 - 5; % correction to have W0 referenced from the 0 longitude(instead of -5E Salih Crater)
dW = 262.7318870466;

% --- Nutation Arguments Sj = Sj0 + Sj_rate * T (degrees) ---
% These correspond to BODY6_NUT_PREC_ANGLES from the SPICE kernel (SSD model part)
% S01 to S18
Sj_constants = [ % [Sj0 (deg), Sj_rate (deg/century)]
    335.844470,   51.7682239;    % S01
    355.351814,  101.6467750;    % S02
      9.369346, 1004.8728024;    % S03
    129.755966, 1223.2050690;    % S04
    219.755966, 1223.2050690;    % S05
    159.835559, 2445.2902118;    % S06
    249.835559, 2445.2902118;    % S07
    117.392885, 3667.0200695;    % S08
    280.169482, 7226.3782354;    % S09
      6.997174, 36506.5422127;   % S10
    196.673251, 15227.2035409;   % S11
    253.848856,  3258.6617087;   % S12
    136.859155,  9266.8742489;   % S13
    144.630256, 12292.3910895;   % S14
      9.821866, 16090.5831593;   % S15
    226.334387, 17383.5986496;   % S16
     93.360491, 18531.0794323;   % S17
     10.9818392,9583937.8056363]; % S18

S_deg = zeros(18, 1);
for j = 1:18
    S_deg(j) = Sj_constants(j,1) + Sj_constants(j,2) * T;
end

% --- Coefficients for nutation series (from SPICE BODY602_NUT_PREC_RA/DEC/PM) ---
% Coefficients Aj for alpha sum(Aj*sin(Sj)) (for S01 to S17)
coeffs_A = [
     0.026616;   % S01
     0.000686;   % S02
    -0.000472;   % S03
    -0.000897;   % S04
     0.002970;   % S05
     0.001127;   % S06
     0.000519;   % S07
     0.000228;   % S08
     0.036804;   % S09
    -0.001107;   % S10
     0.073107;   % S11
    -0.000167;   % S12
     0.000000;   % S13 (from SPICE BODY602_NUT_PREC_RA(21))
     0.000000;   % S14 (from SPICE BODY602_NUT_PREC_RA(22))
    -0.000376;   % S15
     0.000248;   % S16
    -0.000137    % S17
     0.000000    % S18
];

% Coefficients Bj for delta sum(Bj*cos(Sj)) (for S01 to S17)
coeffs_B = [
     0.004398;   % S01
    -0.000264;   % S02
    -0.000185;   % S03
    -0.000093;   % S04
    -0.000068;   % S05
    -0.000236;   % S06
     0.000000;   % S07 (from SPICE BODY602_NUT_PREC_DEC(15))
    -0.000028;   % S08
     0.004141;   % S09
    -0.000124;   % S10
     0.008229;   % S11
     0.000007;   % S12
     0.000000;   % S13 (from SPICE BODY602_NUT_PREC_DEC(21))
     0.000000;   % S14 (from SPICE BODY602_NUT_PREC_DEC(22))
    -0.000039;   % S15
     0.000026;   % S16
    -0.000016    % S17
     0.000000    % S18
];

% Coefficients Cj for W sum(Cj*sin(Sj)) (for S01 to S18)
coeffs_C = [
    -0.026447;   % S01
    -0.000682;   % S02
     0.000469;   % S03
    -0.005118;   % S04
     0.036955;   % S05
    -0.013111;   % S06
     0.014206;   % S07
    -0.006687;   % S08
    -0.036404;   % S09
     0.001082;   % S10
    -0.072604;   % S11
    -0.266358;   % S12
    -0.188429;   % S13
    -0.004710;   % S14
     0.000337;   % S15
    -0.000183;   % S16
    -0.001724;   % S17
    -0.091295    % S18
];

% --- Calculate RA ---
alpha_sum_terms = 0;
for j = 1:18 
    alpha_sum_terms = alpha_sum_terms + coeffs_A(j) * sind(S_deg(j));
end
RA = RA0 + dRA * T + alpha_sum_terms;

% --- Calculate DEC ---
delta_sum_terms = 0;
for j = 1:18
    delta_sum_terms = delta_sum_terms + coeffs_B(j) * cosd(S_deg(j));
end
DEC = DEC0 + dDEC * T + delta_sum_terms;

% --- Calculate W ---
W_sum_terms = 0;
for j = 1:18 
    W_sum_terms = W_sum_terms + coeffs_C(j) * sind(S_deg(j));
end
W = mod(W0 + dW * d + W_sum_terms, 360);

end