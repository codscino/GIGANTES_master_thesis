clc
clear all
close all
%%
% date1 = [2000, 1, 1, 11, 48, 0];
date1 = [2009, 08, 11, 12, 0, 0]; % spring equinox
% date1 = [2017, 05, 23, 12, 0, 0]; % winter solstice
% date1 = [2025 05 06, 12, 0, 0]; % autumn equinox
% date1 = [2032 04 10, 12, 0, 0]; % winter solstice
T = date2mjd2000(date1);

% claudio enceladuys body fixed from paper
[lat_ebf, lon_ebf] = sunSubPointOnEnceladus(T, true);

% enceladus body fixed from radius and orbit
[lat_rad, lon_rad] = sunSubPointOnEnceladusRad(T);


%% spice routine method
encPath = which('fake_enceladus.tpc' );
cspice_furnsh(encPath);

% 3) Compute the geometric sub-solar intercept on a spherical Enceladus
[ spoint, ~, ~ ] = cspice_subslr( ...
    'Intercept/ELLIPSOID', ...                 
    'Enceladus', T*86400, ...
    'IAU_ENCELADUS', ...                      
    'NONE', ...                                 
    'SUN' );                              

% 4) Normalize and convert to planetocentric lat/lon (degrees)
r   = norm( spoint );
lat_spice = asind( spoint(3) / r );                   
lon_spice = atan2d( spoint(2), spoint(1) ) + 5;    % correction to have W0 referenced from the 0 longitude(instead of -5E)      

% 5) Display results
% --- compute errors ---
lat_err = abs( lat_spice - lat_rad );
lon_err = abs( lon_spice - lon_rad);  

% --- print them ---
fprintf( 'error latitude between spice and enceladus radius-orbit rf : %.12f°°\n', lat_err );
fprintf( 'error longitude between spice and enceladus radius-orbit rf : %.12f°\n', lon_err );



%% plot spice vs rad
% define endpoints and step
startDate = [2010,01,01,12,0,0];
endDate   = [2020,01,01,12,0,0];
stepDays  = 1;

% convert to MJD2000
T0 = date2mjd2000(startDate);
T1 = date2mjd2000(endDate);

% build time vector
Ts = T0:stepDays:T1;
N  = numel(Ts);

% preallocate
lat_err_series = zeros(1,N);
lon_err_series = zeros(1,N);

% run the loop
for idx = 1:N
    Tcur = Ts(idx);

    % your radius-orbit method
    [lat_rad, lon_rad] = sunSubPointOnEnceladusRad(Tcur);

    % SPICE intercept
    cspice_furnsh(encPath);
    [spoint, ~, ~] = cspice_subslr( ...
        'Intercept/ELLIPSOID','Enceladus', Tcur*86400, ...
        'IAU_ENCELADUS','NONE','SUN' );

    r = norm(spoint);
    lat_sp = asind(spoint(3)/r);
    lon_sp = atan2d(spoint(2),spoint(1)) + 5;

    % errors
    lat_err_series(idx) = abs(lat_sp - lat_rad);
    lon_err_series(idx) = abs(lon_sp - lon_rad);
end

% plot on logarithmic scale
figure;
semilogy(Ts, lat_err_series, '-o', Ts, lon_err_series, '-s','LineWidth',1.2);
grid on;

% capture handles for labels, legend, title
hXLabel = xlabel('Date');
hYLabel = ylabel('Error (deg)');
hLeg    = legend('Latitude error','Longitude error','Location','SouthEast');
hTitle  = title('Error: SPICE vs Radius–Orbit Frame (1-days step)');

% bump everything to 3× the base axes font
ax       = gca;
baseFS   = ax.FontSize;
hLeg.FontSize    = 3*baseFS;
hTitle.FontSize  = 3*baseFS;
hXLabel.FontSize = 3*baseFS;
hYLabel.FontSize = 3*baseFS;

% also enlarge the tick labels on both axes
ax.FontSize       = 3*baseFS;
ax.XAxis.FontSize = 3*baseFS;
ax.YAxis.FontSize = 3*baseFS;

% 1) Build full year list
dateVecs = zeros(numel(Ts),6);
for i = 1:numel(Ts)
    dateVecs(i,:) = mjd20002date(Ts(i));
end
years = dateVecs(:,1);
uniqueYears = unique(years,'stable');

% 2) Pick every year, beginning with the first
startYear = uniqueYears(1);
yearTicks = startYear : 1 : uniqueYears(end);

% 3) Find the index of the first occurrence of each tick‐year
tickIdx = arrayfun(@(yr) find(years==yr,1,'first'), yearTicks);

% 4) Apply to axes
ax = gca;
ax.XTick = Ts(tickIdx);
ax.XTickLabel = string(yearTicks);
ax.XTickLabelRotation = 45;

%% plot spice vs claudio
% preallocate error series
lat_err_cl = zeros(1,N);
lon_err_cl = zeros(1,N);

for idx = 1:N
    Tcur = Ts(idx);

    % 1) Claudio’s sunSubPointOnEnceladus (paper-based frame)
    [lat_cl, lon_cl] = sunSubPointOnEnceladus(Tcur, true);

    % 2) SPICE intercept on sphere (as before)
    [spoint, ~, ~] = cspice_subslr( ...
        'Intercept/ELLIPSOID','Enceladus', Tcur*86400, ...
        'IAU_ENCELADUS','NONE','SUN' );


    r = norm(spoint);
    lat_sp = asind(spoint(3)/r);
    lon_sp = atan2d(spoint(2),spoint(1)) + 5;

    % 3) Compute errors
    lat_err_cl(idx) = abs(lat_sp - lat_cl);
    lon_err_cl(idx) = abs(lon_sp - lon_cl);
end

% 4) Plot
figure;
semilogy(Ts, lat_err_cl, '-o', Ts, lon_err_cl, '-s','LineWidth',1.2);
grid on;
hXLabel = xlabel('Date');
hYLabel = ylabel('Error (deg)');
hLeg2   = legend('Latitude error','Longitude error','Location','NorthWest');
hTitle2 = title('Error: SPICE vs Claudio’s Orientation (1-day step)');

% Triple the legend, title & axis‐label font sizes
ax = gca;
fs = ax.FontSize;
hLeg2.FontSize   = 3*fs;
hTitle2.FontSize = 3*fs;
hXLabel.FontSize = 3*fs;
hYLabel.FontSize = 3*fs;

% enlarge the tick‐label font to match
ax.FontSize = 3*fs;

% Re‐use year ticks every 5 years
ax.XTick = Ts(tickIdx);
ax.XTickLabel = string(yearTicks);
ax.XTickLabelRotation = 45;