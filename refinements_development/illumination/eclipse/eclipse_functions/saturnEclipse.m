function alpha_saturn = saturnEclipse(beta, d)
    % INPUTS:
    % beta: the sun elevation angle [deg]
    % d: distance between Saturn and Enceladus surface point[ km]
    %
    % OUTPUTS:
    % - alpha_saturn = angular size of saturn [radians]
    
    % equatorial and polar Saturn radii(values from Spice)
    Req = 60268;
    Rpol = 54364;
    beta = deg2rad(beta);

    num = Req * Rpol;
    den1 = (d^2 - Req^2);
    den2 = (Rpol^2 *cos(beta)^2 + Req^2* sin(beta)^2) ;
    alpha_saturn =  atan(num / sqrt( den1 * den2 ));
end