function [] = Plot_Flyby_GT_illuminationv3(Flyby_Data, color, T0, spiceParam, lineStyle, width)
%
% INPUTS:
%   Flyby_Data – struct with fields:
%                .lats     (rad), .longs      (rad),
%                .rp_lat   (deg), .rp_long    (deg),
%   color       – [1×3] RGB for day side
%   T0          – start time in MJD2000 (days)
%   spiceParam  - struct with SPICE parameters
%   lineStyle   - (Optional) Line style ('-', '--', ':', etc.). Default: '-'.
%   width       - (Optional) Line width. Default: 1.5.
%
% Author: adapted from Jose Carlos García Mateas
% Last revision: 2025-07-25 (Merged version)

%—— Set Defaults for Optional Inputs —————————————————
if nargin < 5
    lineStyle = '-'; % Default to a solid line
end
if nargin < 6
    width = 1.5; % Default line width
end

%—— Unpack & Convert Data —————————————————————————————
lats_rad  = Flyby_Data.lats;
longs_rad = Flyby_Data.longs;

% Convert to degrees and wrap longitude to the [0, 360) range
lats  = rad2deg(lats_rad);
longs = mod(rad2deg(longs_rad), 360);

%—— Calculate Illumination at T0 —————————————————
[~, lon_sun_T0] = sunSubPointOnEnceladus(T0, spiceParam); % deg
delta_lon = wrapTo180(longs - lon_sun_T0);
isShade = abs(delta_lon) > 90;

%—— Prepare Data for Discontinuous Plotting ———————————
% Create separate coordinate arrays for day and night segments.
% Replace points of the other type with NaN to create breaks in the line.
day_longs = longs;   day_lats = lats;
day_longs(isShade) = NaN;
day_lats(isShade) = NaN;

night_longs = longs; night_lats = lats;
night_longs(~isShade) = NaN;
night_lats(~isShade) = NaN;

% Find where longitude wraps around to prevent plotting a line across the map.
stop_index = find(diff(longs) > 300, 1);

% If a wrap is found, insert a NaN at that point in all coordinate arrays
% to force a break in the plotted line.
if ~isempty(stop_index)
    day_longs   = [day_longs(1:stop_index),   NaN, day_longs(stop_index+1:end)];
    day_lats    = [day_lats(1:stop_index),    NaN, day_lats(stop_index+1:end)];
    night_longs = [night_longs(1:stop_index), NaN, night_longs(stop_index+1:end)];
    night_lats  = [night_lats(1:stop_index),  NaN, night_lats(stop_index+1:end)];
end

%—— Plot Groundtrack and Special Points —————————————————
hold on;
legend_handles = [];
legend_labels = {};

% Plot the day side track segments. A legend entry is created only if day points exist.
if any(~isShade)
    p_day = plot(day_longs, day_lats, 'LineStyle', lineStyle, 'Color', color, 'LineWidth', width);
    legend_handles(end+1) = p_day;
    legend_labels{end+1} = 'Day side';
end

% Plot the night side track segments. A legend entry is created only if night points exist.
if any(isShade)
    p_night = plot(night_longs, night_lats, 'LineStyle', lineStyle, 'Color', [0 0 0], 'LineWidth', width);
    legend_handles(end+1) = p_night;
    legend_labels{end+1} = 'Night side';
end

% Plot entry, exit, and periapsis markers
p_entry = plot(longs(1), lats(1), 'o', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black', 'MarkerSize', 5);
p_exit = plot(longs(end), lats(end), 'o', 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'MarkerEdgeColor', 'black', 'MarkerSize', 5);
p_peri = plot(Flyby_Data.rp_long, Flyby_Data.rp_lat, 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black', 'MarkerSize', 5);

%—— Finalize Axes & Legend ——————————————————————
xlim([0 360]);

% Combine all handles and labels for a single, accurate legend
all_handles = [legend_handles, p_entry, p_exit, p_peri];
all_labels = [legend_labels, 'Entry', 'Exit', 'Periapsis'];
legend(all_handles, all_labels, 'Location', 'best');

hold off;
end