function beta = saturn_sun_beta_angle(r_sun, r_sat, r_enc, tp)
    % Compute the Sun’s declination from Saturn equatorial in Enceladus
    % Body Fixed reference frame
    %   Input:
    %     tp - time at the pericentre in j2000 days
    %   Output:
    %     beta  - declination angle between Saturn and Sun (deg)
    
    % convert r_sun and r_sat in enceladus body fixed
    u_sat = icrf2enceladus(r_sat, r_enc, tp);
    u_sun = icrf2enceladus(r_sun, r_enc, tp);
    
    % get declinations
    dec_sat = asind(u_sat(3));
    dec_sun = asind(u_sun(3));
    
    beta = abs(dec_sun-dec_sat);

end