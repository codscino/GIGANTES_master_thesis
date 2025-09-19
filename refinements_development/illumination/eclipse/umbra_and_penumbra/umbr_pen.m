function [h_umb, h_pen] = umbr_pen(r_sun, r_planet, d_sp, d_pm)
    % INPUTS:
    % r_sun  = sun radius [km]
    % r_planet = planet radius [km]
    % d_sp = average sun-planet distance [km]
    % d_pm = average planet-moon distance [km]
    %
    % OUTPUTS:
    % h_umb = radius of the umbra [km]
    % h_pen = radius of the penumbra [km]
    
    alpha_umb = rad2deg(asin((r_sun-r_planet)/d_sp));
    alpha_pen = rad2deg(asin((r_sun+r_planet)/d_sp));
    
    x = r_planet/sind(alpha_pen);
    y = r_planet/sind(alpha_umb);
    
    h_umb = (y-d_pm) * tand(alpha_umb);
    h_pen = (x+d_pm) * tand(alpha_pen);

end