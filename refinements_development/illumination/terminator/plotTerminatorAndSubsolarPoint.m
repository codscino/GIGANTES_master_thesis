function plotTerminatorAndSubsolarPoint(T, idCentral, idMoon)
    %Plots the terminator and solar sub-point for a given J2000 date
    %   
    %Input:
    % jd: J2000 date
    % idCentral (optional): id of the central body (default: Saturn, 6).
    % idMoon (optional):id of the moon (default: Enceladus, 1).
    %
    %Warnings:
    % 1. The proper SPICE kernels must be loaded before launching the script
    % to make sunSubPointOnEnceladus function work.
    %
    % 2. At the moment the script works just with Saturn and Enceladus because
    % sunSubPointOnEnceladus has not been genralized


    % Set default values for idCentral and idMoon
    if nargin < 3 || isempty(idMoon)
        idMoon = 1;    % Default Enceladus
    end
    if nargin < 2 || isempty(idCentral)
        idCentral = 6; % Default Saturn
    end


    % Compute the solar sub-point
    [lat_sun, lon_sun] = sunSubPointOnEnceladus(T);
    % Ensure longitude is [0, 360) for consistent plotting
    lon_sun = mod(lon_sun + 360, 360);

    % Draw the base Mercator texture
    plotTextureLatLong(idMoon, idCentral);
    hold on; % Keep the plot active for adding more elements

    %% Compute and plot the terminator curve
    npts = 361; % resolution for the curve
    % Define longitudes for plotting the terminator curve
    lon_plot = linspace(0, 360, npts);

    lha = deg2rad(lon_plot - lon_sun); % local hour angle [rad]
    
    % Formula for terminator latitude
    lat_formula = atan( -cos(lha) ./ tan(deg2rad(lat_sun)) ); % terminator latitude [rad]
    lat = rad2deg(lat_formula); % convert back to degrees

    %% Shade the nightside
    % Nightside is the region towards the pole opposite to the sub-solar latitude.
    if lat_sun >= 0
        edgeLat = -90;  % extend patch down to South Pole
    else
        edgeLat = +90;  % extend patch up   to North Pole
    end

    px = [lon_plot, fliplr(lon_plot)];
    py = [lat, edgeLat*ones(1, npts)];
    patch(px, py, 'k', 'FaceAlpha', 0.4, 'EdgeColor', 'none');

    %% Draw the terminator line
    hTerm = plot(lon_plot, lat, 'r-', 'LineWidth', 2);

    %% Highlight the sub-solar point
    hSun = plot(lon_sun, lat_sun, 'yo', 'MarkerFaceColor', ...
        'y', 'MarkerSize', 8);

    %% Extras
    xlim([0 360])

    legend([hTerm, hSun], {'Terminator','Sub-solar point'}, ...
        'Location','best')
    title( sprintf(['Terminator Shadow (\\lambda_{sun}=%.1f°, ' ...
        '\\phi_{sun}=%.1f°)'], lon_sun, lat_sun) )

    hold off; 
end