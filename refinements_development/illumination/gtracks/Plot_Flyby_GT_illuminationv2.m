function [] = Plot_Flyby_GT_illuminationv2(Flyby_Data, color, T0, spiceParam)
% Plot_Flyby_GT_illumination  Plots a flyby groundtrack with day/night shading.
%
% USAGE:
%   Plot_Flyby_GT_illumination(Flyby_Data, color, T0)
%
% INPUTS:
%   Flyby_Data – struct with fields:
%                .lats     (rad), .longs      (rad),
%                .rp_lat   (deg), .rp_long    (deg),
%   color       – [1×3] RGB for day side
%   T0          – start time in MJD2000 (days)
%
% Author: adapted from Jose Carlos García Mateas
% Last revision: 2025-05-26

%—— unpack & convert —————————————————————————————
lats_rad  = Flyby_Data.lats;
longs_rad = Flyby_Data.longs;
N         = numel(lats_rad);

% time vector over 1 hour
ts        = linspace(T0, T0 + 1/24, N);

% convert to degrees and wrap [0,360)
lats      = rad2deg(lats_rad);
longs     = mod(rad2deg(longs_rad),360);


%—— fast illumination test via terminator Great Circle —————————————————
[~, lon_sun0] = sunSubPointOnEnceladus(T0, true, spiceParam);       % deg
[~, lon_sun1] = sunSubPointOnEnceladus(T0+1/24, true, spiceParam);       % deg

dPsi_dt = (lon_sun1 - lon_sun0) / (1/24);   % deg per day
isDayFast = false(1, N);
for i = 1:N
    lon = interp1([T0, T0+1/24], [lon_sun0, lon_sun1], ts(i));
    dLon = wrapTo180(longs(i) - lon);
    isDayFast(i) = abs(dLon) <= 90;
end
isShade = ~isDayFast;

%—— find terminator crossing indices —————————————————————————
crossIdx = find(diff(isDayFast) ~= 0);

%—— refine crossing times with Newton-Raphson —————————————————————
dLon_dt = diff(longs) ./ diff(ts); % deg per day (vector length N-1)
for idx = crossIdx
    t_low = ts(idx);
    t_high= ts(idx+1);
    % determine initial sign & function
    sign0 = double(isDayFast(idx)) * 2 - 1; % +1 for day→night, -1 for night→day reversal
    % initial guess
    t = (t_low + t_high)/2;
    % Newton iterations
    for iter = 1:5
        % interpolate lon & psi
        lon_t = interp1(ts, longs, t);
        psi_t = interp1([T0, T0+1/24], [lon_sun0, lon_sun1], t);
        delta = wrapTo180(lon_t - psi_t);
        F = delta - sign0*90;
        % derivative dF/dt = dlon/dt - dpsi/dt
        % find nearest index for dlon_dt
        k = idx;
        dlon = dLon_dt(k);
        dF = dlon - dPsi_dt;
        t = t - F / dF;
    end
    % evaluate exact shading at crossing
    lat_t = interp1(ts, lats, t);
    lon_t = interp1(ts, longs, t);
    shadeCross = isGTinShade(t, lat_t, lon_t);
    % assign shading flags around crossing
    isShade(idx)   = shadeCross;
    isShade(idx+1) = shadeCross;
end

%—— detect wrap-around in longitude ———————————————————
stop_index = find(diff(longs) > 300, 1);

%—— plot groundtrack with shading —————————————————————
hold on;
if ~isempty(stop_index) && stop_index > 1
    L1 = 1:stop_index;   L2 = stop_index+1 : N;
    plot(longs(L1 & ~isShade(L1)), lats(L1 & ~isShade(L1)), '.', 'Color', color, 'MarkerSize',8);
    plot(longs(L1 & isShade(L1)),  lats(L1 & isShade(L1)),  '.', 'Color',[0 0 0], 'MarkerSize',8);
    plot(longs(L2 & ~isShade(L2)), lats(L2 & ~isShade(L2)), '.', 'Color', color, 'MarkerSize',8);
    plot(longs(L2 & isShade(L2)),  lats(L2 & isShade(L2)),  '.', 'Color',[0 0 0], 'MarkerSize',8);
else
    plot(longs(~isShade), lats(~isShade), '.', 'Color', color, 'MarkerSize',8);
    plot(longs(isShade),  lats(isShade),  '.', 'Color',[0 0 0], 'MarkerSize',8);
end

%—— entry, exit, periapsis points ————————————————
plot(longs(1),    lats(1),    'o', 'MarkerFaceColor','blue',  'MarkerEdgeColor','black','MarkerSize',5);
plot(longs(end),  lats(end),  'o', 'MarkerFaceColor',[0.9290 0.6940 0.1250], 'MarkerEdgeColor','black','MarkerSize',5);
plot(Flyby_Data.rp_long, Flyby_Data.rp_lat, 'o', 'MarkerFaceColor','red', 'MarkerEdgeColor','black','MarkerSize',5);

%—— finalize axes & legend ——————————————————————
xlim([0 360]);
legend('Day side','Night side','Entry','Exit','Periapsis','Location','best');
hold off;
end