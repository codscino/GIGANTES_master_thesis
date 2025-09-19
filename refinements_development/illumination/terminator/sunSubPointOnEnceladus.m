function [lat_sun, lon_sun] = sunSubPointOnEnceladus(T, spiceParam)
    % INPUTs:
    %   T in MJD2000
    %
    % OUTPUTs:
    %   lat_sun, lon_sun -> lat, lon of the subsolar point on Enceladus

    
    if nargin < 2 || isempty(spiceParam)
        spiceParam.frame    = 'J2000';
        spiceParam.abcorr   = 'NONE';
        spiceParam.observer = '0';
    end
    
    r_sun = EphSS_car_spice2(10,T,true, spiceParam);
    r_enc = EphSS_car_spice2(602,T,true, spiceParam);
    
    % sun position versor in ENceldus body fixed reference frame
    u = icrf2enceladus(r_sun, r_enc, T);

    % from  Enceladus body fixed RF to latitude,longitude 
    % assuming a perfect spherical Enceladus(2% oblateness in reality)
    lat_sun = asin(u(3));
    lon_sun = atan2(u(2), u(1));

    lat_sun = rad2deg(lat_sun);
    lon_sun = rad2deg(lon_sun);

end
