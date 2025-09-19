function isShade = isGTinShade(T, lat_gt, lon_gt, spiceParam)
% isGTinShade  True if a given lat/lon on Enceladus is in night.
%
% USAGE:
%   isShade = isGTinShade(T, lat_gt, lon_gt)
%   isShade = isGTinShade(T, lat_gt, lon_gt, isSpice, kernels)
%
% INPUTS:
%   T        - time in MJD2000
%   lat_gt   - target latitude (deg), scalar or vector
%   lon_gt   - target longitude (deg), same size as lat_gt
%   isSpice  - (opt) use SPICE? default = true
%   kernels  - (opt) SPICE kernels cell array; defaults to {'sat441.bsp','naif0012.tls'}
%
% OUTPUT:
%   isShade  - logical array, true if in night (terminator‐based)

    if nargin < 4
        spiceParam.frame    = 'J2000';
        spiceParam.abcorr   = 'NONE';
        spiceParam.observer = '0';
    end

    [lat_sun, lon_sun] = sunSubPointOnEnceladus(T, true, spiceParam);

    % bring longitudes into (0,360)
    lon_sun = mod(lon_sun,360);
    lon_gt  = mod(lon_gt,360);

    % compute local hour angle of each gt point relative to sub-solar
    lha = deg2rad(lon_gt - lon_sun);

    % terminator latitude at each lon_gt (deg)
    lat_term = atan( -cos(lha) ./ tan(deg2rad(lat_sun)) );
    lat_term = rad2deg(lat_term);

    % night is the side away from the subsolar latitude
    if lat_sun >= 0
        % Sun to north: night is below terminator
        isShade = (lat_gt < lat_term);
    else
        % Sun to south: night is above terminator
        isShade = (lat_gt > lat_term);
    end
end