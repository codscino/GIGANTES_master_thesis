function alpha = getAngSize(r_observer, r_target, id_target)
    % function to find the angular size of a target spherical body seen 
    % from an observer body. 
    %
    % The distance between the 2 bodies is considere much bigger than the
    % radius of the observer, so the observation is approximated from the
    % centre of the observer body.
    % The target body is approximated as a sphere
    % 
    % INPUTS:
    % - id_observer -> NAIF_id integer of the observer body(ex. 602 for Enceladus)
    % - id_target -> NAIF_id integer of the target body(ex. 10 for Sun)
    % - T -> epoch in mjd2000 days
    %
    % OUTPUTS:
    % - alpha = angular size in radians
    
    d = norm(r_target-r_observer);
    
    if id_target == 10 % default value for the Sun
        R_target = 696000;
    else
        [radii_target] = cspice_bodvrd(num2str(id_target), 'RADII', 3);
        R_target = mean(radii_target); % target body approximated as a sphere
    end
    
    alpha = asin(R_target/d);
end