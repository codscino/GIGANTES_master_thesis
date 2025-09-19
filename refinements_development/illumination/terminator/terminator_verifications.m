clc;
clear all;
close all;

% a script to compute multiple verifications on terminator plot logic
% developed in Claudio Ferrara master thesis

% inputs
idCentral = 6;
idMoon    = 1;

kernels = {'de440.bsp', 'sat441.bsp', 'naif0012.tls'};
loadSpiceKernels(kernels);

% https://www.planetary.org/articles/06031044-oppositions-conjunctions-rpx 
spring_equinox = [2009, 08, 11, 12, 0, 0];
summer_solstice = [2017, 05, 23, 12, 0, 0];
autumn_equinox = [2025 05 06, 12, 0, 0];
winter_solstice = [2032 04 10, 12, 0, 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Rectangular terminator on spring equinox %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% qualitative proof that at equinox the solar sub point latitude is zero
% and the terminator becomes a staright vertical line

se = date2mjd2000(spring_equinox); % convert spring equinox in julian days
plotTerminatorAndSubsolarPoint(se);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Summer solstice terminator %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% qualitative proof that at summer solstice the solar sub-point latitude 
% is close to 26.72 deg, the inclination of Saturn(Enceladus orbit coplanar 
%to Saturn)

ss = date2mjd2000(summer_solstice);
plotTerminatorAndSubsolarPoint(ss);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Solar sub point verifications %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Comparison of the sub-solar point latitude and longitude using:
% A) Original Claudio script
% B) SPICE with a modified spherical Enceladus kernel
% C) SPICE with original ellipsoidal ENceladus kernel


T = date2mjd2000(summer_solstice);
% use SPICE time to be more accurate
% T_spice = UTC2ET(T); 

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% A. Original script %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate solar sub-point using script developed by Claudio
[lat_sun_claudio, lon_sun_claudio] = sunSubPointOnEnceladus(T);
lon_sun_claudio = mod(lon_sun_claudio + 360, 360);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% B. Spherical Enecleadus %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loadSpiceKernels({'fake_enceladus.tpc'}); % modified kernel with all
                                          % radii = 252 km

% Compute the geometric sub-solar intercept on a spherical Enceladus
[ spoint_spherical, ~, ~ ] = cspice_subslr( ...
    'NEAR POINT/ELLIPSOID', ...                 
    'Enceladus', T*86400, ...
    'IAU_ENCELADUS', ...                     
    'NONE', ...                               
    'SUN' );    

% Normalize and convert to planetocentric lat/lon (degrees)
r_spherical   = norm( spoint_spherical );
lat_sun_spherEnc = asind( spoint_spherical(3) / r_spherical );                          
lon_sun_spherEnc = atan2d( spoint_spherical(2), spoint_spherical(1) );        
lon_sun_spherEnc = mod(lon_sun_spherEnc + 360, 360);

unloadSpiceKernels({'fake_enceladus.tpc'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% C. Ellipsoidal Enecleadus %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loadSpiceKernels({'enceladus.tpc'}); % real kernel with different radii

% Compute the geometric sub-solar intercept on a ellipsoidal Enceladus
[ spoint_ellipsoidal, ~, ~] = cspice_subslr( ...
    'NEAR POINT/ELLIPSOID', ...                 
    'Enceladus', T*86400, ...
    'IAU_ENCELADUS', ...                     
    'NONE', ...                               
    'SUN' );

unloadSpiceKernels({'enceladus.tpc'})

% Normalize and convert to planetodetic lat/lon (degrees)
re =  256.14;
rp =  248.68 ;
f = (re - rp) / re;

% use spice routine for ellipsoidal latitute and logitudes (planetodetic)
[lon_rad_geo, lat_rad_geo] = cspice_recgeo(spoint_ellipsoidal, re, f);
lat_sun_ellipsEnc = rad2deg(lat_rad_geo);                         
lon_sun_ellipsEnc = rad2deg(lon_rad_geo);       
lon_sun_ellipsEnc = mod(lon_sun_ellipsEnc + 360, 360);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% -- Solar sub-point comaprison plot -- %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotTextureLatLong(idMoon, idCentral);
hold on; 

% Plot the three sub-solar points
hClaudio = plot(lon_sun_claudio, lat_sun_claudio, ...
    'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
hSphere  = plot(lon_sun_spherEnc, lat_sun_spherEnc, ...
    'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 6); % Square marker for spherical
hEllipso = plot(lon_sun_ellipsEnc, lat_sun_ellipsEnc, ...
    'gd', 'MarkerFaceColor', 'g', 'MarkerSize', 6); % Diamond marker for ellipsoidal

% Add Legend and Title
legend([hClaudio, hSphere, hEllipso], ...
    {'Claudio''s Method', 'Spherical Enceladus (fake kernel)', 'Ellipsoidal Enceladus (real kernel)'}, ...
    'Location', 'best');

% Format the date for the title
title_date_vec = mjd20002date(T);
title_date_str = sprintf('%04d-%02d-%02d %02d:%02d:%02.0f', ...
                         title_date_vec(1), title_date_vec(2), title_date_vec(3), ...
                         title_date_vec(4), title_date_vec(5), title_date_vec(6));

title(sprintf('Sub-solar Point Comparison on Enceladus (Epoch: %s)', title_date_str));

xlim([0 360]);

hold off; 

% Print Error Information
fprintf('\n--- Sub-solar Point Error Comparison (Epoch: %s) ---\n', title_date_str);
fprintf('Claudio vs. Spherical Enceladus:\n');
fprintf('  Latitude Error:  %.4f degrees\n', abs(lat_sun_claudio - lat_sun_spherEnc));
fprintf('  Longitude Error: %.4f degrees\n', abs(lon_sun_claudio - lon_sun_spherEnc));
fprintf('\nClaudio vs. Ellipsoidal Enceladus:\n');
fprintf('  Latitude Error:  %.4f degrees\n', abs(lat_sun_claudio - lat_sun_ellipsEnc));
fprintf('  Longitude Error: %.4f degrees\n', abs(lon_sun_claudio - lon_sun_ellipsEnc));
fprintf('-----------------------------------------------------\n');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Umbral and Penumbral terminator comparison using SPICE %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Qualitative verification to show that the the sun rays could be assumed 
% parallels(infinite distance).
% The plot infact shows there is a negligible differnce between the umbral
% and penumral terminator using the SPICE routine cspice_edterm


loadSpiceKernels({'pck00010.tpc'});
% Define parameters for the terminator calculation
target = 'ENCELADUS';
illuminator = 'SUN';
abcorr = 'NONE'; % Aberration correction
observer = 'ENCELADUS'; % Observer is the target itself for terminator calculation
fixref = 'IAU_ENCELADUS'; % Body-fixed frame for Enceladus
et = cspice_str2et('2025-08-26T12:00:00'); % Example epoch
npts = 361; % Number of points to calculate on the terminator line

% Calculate the UMBRAL terminator points
[trgepc_umbral, obspos_umbral, trmpts_umbral] = cspice_edterm('UMBRAL', illuminator, ...
                                 target, et, fixref, abcorr, observer, npts);

% Convert Cartesian coordinates to latitude and longitude for UMBRAL
[radius_umbral, longitude_umbral_rad, latitude_umbral_rad] = cspice_reclat(trmpts_umbral);
longitude_deg_umbral = mod(rad2deg(longitude_umbral_rad) + 360, 360);
latitude_deg_umbral = rad2deg(latitude_umbral_rad);

% Calculate the PENUMBRAL terminator points
[trgepc_penumbral, obspos_penumbral, trmpts_penumbral] = cspice_edterm('PENUMBRAL', illuminator, ...
                                 target, et, fixref, abcorr, observer, npts);

% Convert Cartesian coordinates to latitude and longitude for PENUMBRAL
[radius_penumbral, longitude_penumbral_rad, latitude_penumbral_rad] = cspice_reclat(trmpts_penumbral);
longitude_deg_penumbral = mod(rad2deg(longitude_penumbral_rad) + 360, 360);
latitude_deg_penumbral = rad2deg(latitude_penumbral_rad);


%%% plot Mercator projection with updated styles %%%
plotTextureLatLong(idMoon, idCentral);
hold on

% Plot Penumbral terminator in a lighter red solid (on top, simulating semi-transparency)
hPenumbral = plot(longitude_deg_penumbral, latitude_deg_penumbral, 'Color', [1 0.5 0.5], 'LineWidth', 2);

% Plot Umbral terminator in black dashed (below)
hUmbral = plot(longitude_deg_umbral, latitude_deg_umbral, 'k--', 'LineWidth', 2);

xlim([0 360])


% Set font size for x and y labels
hXLabel = xlabel('Longitude (deg)');
hYLabel = ylabel('Latitude (deg)');

% Set font size for the legend
hLegend = legend([hUmbral, hPenumbral], {'Umbral Terminator (Dashed Black)', 'Penumbral Terminator (Lighter Red)'},'Location','best');

grid on;
hold off


