function [] = Plot_Flyby_GT_illumination(Flyby_Data, color, T0, timeStep)
% Plot_Flyby_GT_illumination  Plots a flyby groundtrack with day/night shading.
%
% USAGE:
%   Plot_Flyby_GT_illumination(Flyby_Data, color, T0, timeStep)
%
% INPUTS:
%   Flyby_Data – struct with fields:
%                .lats     (rad), .longs      (rad),
%                .rp_lat   (deg), .rp_long    (deg),
%   color       – [1×3] RGB for day side
%   T0          – start time in MJD2000 (days)
%   timeStep    – time step in minutes for shading calculation (default: 1)
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

%—— handle input default ——————————————————————————
if nargin < 4
    timeStep = 1; % 1 minute
end

%—— compute shade flag once per timeStep —————————————————
isShade = false(1,N);
% number of intervals in one hour
totalMin = 60;
numSteps = ceil(totalMin / timeStep);
for k = 1 : numSteps
    t0      = T0 + (k-1)*timeStep/1440;
    t1      = T0 + min(k*timeStep, totalMin)/1440;
    idx0    = find(ts >= t0 & ts < t1, 1, 'first');
    if isempty(idx0), continue; end
    shadeFlag = isGTinShade(ts(idx0), lats(idx0), longs(idx0));
    % apply to all samples in this interval
    intervalIdx = (ts >= t0) & (ts < t1);
    isShade(intervalIdx) = shadeFlag;
end

%—— detect wrap-around in longitude ———————————————————
stop_index = find(diff(longs) > 300, 1);

%—— plot groundtrack with shading —————————————————————
if ~isempty(stop_index) && stop_index > 1
    L1 = 1:stop_index;   L2 = stop_index+1 : N;
    % segment 1
    plot(longs(L1 & ~isShade(L1)), lats(L1 & ~isShade(L1)), '.', 'Color', color, 'MarkerSize',8); hold on;
    plot(longs(L1 & isShade(L1)),  lats(L1 & isShade(L1)),  '.', 'Color',[0 0 0], 'MarkerSize',8);
    % segment 2
    plot(longs(L2 & ~isShade(L2)), lats(L2 & ~isShade(L2)), '.', 'Color', color, 'MarkerSize',8);
    plot(longs(L2 & isShade(L2)),  lats(L2 & isShade(L2)),  '.', 'Color',[0 0 0], 'MarkerSize',8);
else
    plot(longs(~isShade), lats(~isShade), '.', 'Color', color, 'MarkerSize',8); hold on;
    plot(longs(isShade),  lats(isShade),  '.', 'Color',[0 0 0], 'MarkerSize',8);
end

%—— entry, exit, periapsis points ————————————————
plot(longs(1),    lats(1),    'o', 'MarkerFaceColor','blue',  'MarkerEdgeColor','black','MarkerSize',5);
plot(longs(end),  lats(end),  'o', 'MarkerFaceColor',[0.9290 0.6940 0.1250], 'MarkerEdgeColor','black','MarkerSize',5);
plot(Flyby_Data.rp_long, Flyby_Data.rp_lat, 'o', 'MarkerFaceColor','red', 'MarkerEdgeColor','black','MarkerSize',5);

%—— finalize axes & legend ——————————————————————
xlim([0 360])
legend('Day side','Night side','Entry','Exit','Periapsis','Location','best')
title(sprintf('Flyby Groundtrack (%d pts over 1 h, \Delta=%g min)', N, timeStep))
hold off
end
