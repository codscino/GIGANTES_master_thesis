function nu_deg = getEnceladusTrueAnomaly(T)
% getEnceladusTrueAnomaly  Compute Enceladus’s true anomaly [deg] via CSPICE
%
% INPUTS:
%   T       = epoch in MJD2000 (days since J2000.0)
%   kernels = cell array of exactly these three SPICE kernels:
%             {
%               'naif0012.tls',  ...  % leap seconds
%               'pck00010.tpc',  ...  % planetary constants (Saturn GM, shapes, etc.)
%               'sat441.bsp'     ...  % Saturn‐system SPK (Enceladus→Saturn)
%             }
%
% OUTPUT:
%   nu_deg = true anomaly (degrees) of Enceladus about Saturn at time T
%
% STEPS:
%  1) Locate each kernel file via which(), then load (cspice_furnsh).
%  2) Convert T to ephemeris seconds (ET).
%  3) Use cspice_spkezr to get Enceladus’s state [r; v] w.r.t. Saturn.
%  4) Call cspice_bodvrd('SATURN','GM',1) to get Saturn’s GM.
%  5) Call cspice_oscelt on [r; v], ET, and GM to get [a, e, i, Ω, ω, M].
%  6) Compute true anomaly via nu = M2theta(M, e) and convert to degrees.
%  7) Unload kernels.
        

    % 3) Convert T (days since J2000) to ephemeris time (seconds past J2000)
    et = T * 86400;

    % 4) Get Enceladus’s state relative to Saturn (ID=6) in ICRF93, no light‐time correction
    [state602, ~] = cspice_spkezr('602', et, 'J2000', 'NONE', '6');
    % state602 is 6×1: [r_x; r_y; r_z; v_x; v_y; v_z]

    % 5) Fetch Saturn’s GM from the loaded PCK
    gmSat = getAstroConstants('Saturn','Mu');

    % 6) Recover osculating Kepler elements: [a, e, i, Ω, ω, M, ...]
    elts = cspice_oscelt(state602, et, gmSat);
    % elts(1)=semimajor axis [km], elts(2)=eccentricity, elts(6)=mean anomaly [rad]
    e  = elts(2);
    M0 = elts(6);

    % 7) Compute true anomaly using M2theta (assumed on the path):
    nu = M2theta(M0, e);      % returns ν in radians
    nu_deg = rad2deg(nu);     % convert to degrees
end