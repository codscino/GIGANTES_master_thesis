function angle_sep = angle_sep(lat_enceladus, lon_enceladus, T, r_enceladus_icrf, r_target_icrf, r_sun_icrf, R_enc)
%
% This function calculates the angular separation between
% a target body(ex. Saturn) and the Sun as viewed from a specific location
% on the surface of Enceladus at a given time T
%
%  Author        : Claudio Ferrara
%
%  Inputs:
%    lat_enceladus    - latitude of the site on Enceladus  [rad]
%    lon_enceladus    - Longitude of the site on Enceladus          [rad]
%    T                - J2000 ET in days, time of the observation
%    r_enceladus_icrf - Position vector of Enceladus' center in ICRF [km] (1x3 row vector)
%    r_target_icrf    - Position vector of the target body's center in ICRF [km] (1x3 row vector)
%    r_sun_icrf       - Position vector of the Sun's center in ICRF  [km] (1x3 row vector)
%    R_enc            - Mean radius of Enceladus                     [km]
%
%  Outputs:
%    angle_sep        - Angular separation between target body and Sun [rad]
%

    % --- Step 1: Calculate the observer's position in ICRF (once) ---
    
    % Convert observer's lat/lon to a cartesian unit vector in the
    % Enceladus Body-Fixed (EBF) frame, as a row vector.
    u_ebf = [cos(lat_enceladus) * cos(lon_enceladus), ...
             cos(lat_enceladus) * sin(lon_enceladus), ...
             sin(lat_enceladus)];
    
    % Get the absolute position of the observer in ICRF 
    r_obs_icrf = enceladus2icrf(u_ebf, r_enceladus_icrf, R_enc, T);
  
    % --- Step 2: Find the topocentric vectors from the observer to each body ---
    rho_target_icrf = r_target_icrf - r_obs_icrf; % from observer to target body
    rho_sun_icrf    = r_sun_icrf    - r_obs_icrf; % from observer to Sun
    % rho_target_icrf = r_target_icrf - r_enceladus_icrf; % from observer to target body
    % rho_sun_icrf    = r_sun_icrf    - r_enceladus_icrf; % from observer to Sun
    
    % --- Step 3: Calculate the angular separation using the dot product ---
    cos_angle = dot(rho_target_icrf, rho_sun_icrf) / (norm(rho_target_icrf) * norm(rho_sun_icrf));
    
    % Clamp the value to [-1, 1] to prevent floating point errors from
    % causing acos to return a complex number 
    angle_sep = acos(max(-1.0, min(1.0, cos_angle)));

end