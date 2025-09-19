function [lat, lon] = sunSubPointOnEnceladusRad(T, kernels)
    % sunSubPointOnEnceladusSpice  Sub-solar lat/lon on spherical Enceladus
    %  
    % INPUT:
    %   T        = time in MJD2000
    %   kernels  = cell array, e.g. {'sat441.bsp','naif0012.tls'}
    % OUTPUT (degrees):
    %   lat, lon = planetocentric latitude & longitude of the sun
    %
    % FRAME: 
    %   X_e → toward Saturn
    %   Z_e → orbit normal (pos × vel)
    %   Y_e = Z_e × X_e
    
    if nargin < 2
        SPKfile           = 'sat441.bsp';
        leapsecondsKernel = 'naif0012.tls';
        kernels           = {SPKfile, leapsecondsKernel};
    end


    % 1) Get Sun and Enceladus states (position & velocity) in SSB/ICRF
    [r_sun, ~]      = EphSS_car_spice(10,  T, true, kernels);  % Sun wrt SSB
    [r_enc, v_enc]  = EphSS_car_spice(602, T, true, kernels);  % Enceladus wrt SSB
    [r_sat, v_sat]      = EphSS_car_spice(6,   T, true, kernels);  % Saturn wrt SSB

    % 2) Sub-solar direction in ICRF
    u_se = (r_sun - r_enc);
    u_se = u_se / norm(u_se);
    u_se = u_se';

    % 3) Build Enceladus‐fixed axes in ICRF

    % 3a) X_e: from Enc toward Saturn
    r_es = r_sat - r_enc;
    x_e = (r_es / norm(r_es))';

    % 3b) Z_e: orbit normal (pos × vel)
    v_es = v_sat - v_enc; 
    z_e = cross( r_es, v_es);
    z_e = (z_e / norm(z_e))';

    % 3c) Y_e: completes RH triad
    y_e = cross( z_e, x_e );
    y_e = y_e / norm(y_e);

    % 4) Rotation matrix ICRF→Enc‐fixed: columns are [x_e y_e z_e]
    C = [ x_e, y_e, z_e ];

    % 5) Express sun‐direction in body frame
    v_bf = C' * u_se;

    % 6) Spherical lat/lon on unit sphere
    lat = asind( v_bf(3) );        % planetocentric latitude
    lon = atan2d( v_bf(2), v_bf(1) ); % planetocentric longitude

end
