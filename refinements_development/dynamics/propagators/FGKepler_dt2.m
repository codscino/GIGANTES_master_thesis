function [rf, vf] = FGKepler_dt2(r0,v0,dt,mu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FGKepler_dt2 - Compute final state vectors (rf and vf) after a given dt
%
%   Inputs:
%       r0  - Initial position vector (km)
%       v0  - Initial velocity vector (km/s)
%       dt  - Time interval over which the state is propagated (seconds)
%       mu  - Gravitational parameter of the central body (km^3/s^2)
%
%   Outputs:
%       rf  - Final position vector after time dt (km)
%       vf  - Final velocity vector after time dt (km/s)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % be careful this function crashes if a<0
    % check this condition before launching it
    
    R0 = norm(r0);
    V0 = norm(v0);
    a = mu / ((2*mu)/R0 - V0^2);

    
    n = sqrt(mu/a^3);
    dM = n*dt;
    sigma0 = dot(r0,v0) / sqrt(mu);

    % find dE with fzero
    fun = @(dE) -dM + dE - (1- (R0/a))*sin(dE) - (sigma0/sqrt(a))*(cos(dE)-1);
    
    try
        x = fzero(fun, dM);
    catch ME
        fprintf('Error in fzero with dM = %f, fun(dM) = %f, a = %f\n', dM, fun(dM), a);
        rethrow(ME);
    end

    dE = x;

    F = 1 - (a/R0)*(1-cos(dE));
    G = dt + sqrt(a^3/mu)*(sin(dE)-dE);

    rf = F*r0 + G*v0;
    
    % II part
    Rf = norm(rf);
    dG = 1 - (a/Rf)*(1-cos(dE));
    vf = 1/G * (dG*rf - r0);

end