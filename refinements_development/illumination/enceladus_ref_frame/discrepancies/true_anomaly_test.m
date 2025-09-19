clc
clear all
close all

%% kernels
kernels = { ...
            'naif0012.tls', ...
            'pck00010.tpc', ...
            'sat441.bsp' ...
        };
  

% Resolve each kernel’s full path via which()
lskFile = which(string(kernels{1}));
pckFile = which(string(kernels{2}));
spkFile = which(string(kernels{3}));

% 2) Furnish (load) the three kernels
cspice_furnsh({lskFile, pckFile, spkFile});
encPath = which('fake_enceladus.tpc' );
cspice_furnsh(encPath);

%% plot spice vs rad
% define endpoints and step
startDate = [2015,01,1,13,0,0];
endDate   = [2015,2,30,12,0,0];
stepDays  = 0.02;

% convert to MJD2000
T0 = date2mjd2000(startDate);
T1 = date2mjd2000(endDate);

% build time vector
Ts = T0:stepDays:T1;
% Preallocate arrays
N = numel(Ts);
lat_err_series = zeros(1, N);
lon_err_series = zeros(1, N);
nu_series      = zeros(1, N);

for idx = 1:N
    Tcur = Ts(idx);

    % 1) Your radial‐frame subsolar method
    [lat_rad, lon_rad] = sunSubPointOnEnceladusRad(Tcur);

    % 2) SPICE intercept (IAU_ENCELADUS frame)
    [spoint, ~, ~] = cspice_subslr( ...
        'Intercept/ELLIPSOID', 'Enceladus', Tcur*86400, ...
        'IAU_ENCELADUS', 'NONE', 'SUN' );

    r = norm(spoint);
    lat_sp = asind(spoint(3)/r);
    lon_sp = atan2d(spoint(2), spoint(1)) + 5;  % +5° to align crater‐Salih

    % 3) Compute errors
    lat_err_series(idx) = abs(lat_sp - lat_rad);
    lon_err_series(idx) = abs(lon_sp - lon_rad);

    % 4) Compute true anomaly (degrees)
    nu_series(idx) = getEnceladusTrueAnomaly(Tcur);
end

% Filter out longitude‐error values > 2 by setting them to NaN
lon_err_filtered = lon_err_series;
lon_err_filtered(lon_err_filtered > 2) = NaN;

figure;
% original plot with true‐anomaly
% semilogy(Ts, lat_err_series, '-o', ...
%          Ts, lon_err_filtered, '-s', ...
%          Ts, nu_series, '-^', 'LineWidth', 1.2);

semilogy(Ts, lat_err_series, '-o', ...
         Ts, lon_err_filtered, '-s', 'LineWidth', 1.2);

hold on;

% get base font‐size for doubling
ax     = gca;
baseFS = ax.FontSize;

% true anomaly blue lines (thicker + doubled‐size labels in one call)
h1 = yline( 90,  'b-', 'LineWidth', 2, 'Label', '90 deg',  ...
           'FontSize', 2*baseFS, ...
           'LabelHorizontalAlignment','right', ...
           'LabelVerticalAlignment','bottom');

h2 = yline(180,  'b-', 'LineWidth', 2, 'Label', '180 deg', ...
           'FontSize', 2*baseFS, ...
           'LabelHorizontalAlignment','right', ...
           'LabelVerticalAlignment','bottom');

h3 = yline(270,  'b-', 'LineWidth', 2, 'Label', '270 deg', ...
           'FontSize', 2*baseFS, ...
           'LabelHorizontalAlignment','right', ...
           'LabelVerticalAlignment','bottom');

h4 = yline(360,  'b-', 'LineWidth', 2, 'Label', '360 deg', ...
           'FontSize', 2*baseFS, ...
           'LabelHorizontalAlignment','right', ...
           'LabelVerticalAlignment','top');

dk = [0 0.5 0];

% prepare storage for triangle‐points
redX   = [];  redY   = [];
greenX = [];  greenY = [];

% Find local minima and maxima …
for i = 2:(N-1)
    val  = lon_err_filtered(i);
    if isnan(val), continue; end
    prev = lon_err_filtered(i-1);
    next = lon_err_filtered(i+1);
    if ~isnan(prev) && ~isnan(next)
        % local minimum → red down‐triangle
        if val < prev && val < next
            redX(end+1)   = Ts(i);
            redY(end+1)   = nu_series(i);
            semilogy(Ts(i), nu_series(i), 'v', ...
                     'MarkerSize',     8, ...
                     'MarkerEdgeColor','r', ...
                     'MarkerFaceColor','r');
        end
        % local maximum → green up‐triangle
        if val > prev && val > next
            greenX(end+1) = Ts(i);
            greenY(end+1) = nu_series(i);
            semilogy(Ts(i), nu_series(i), '^', ...
                     'MarkerSize',     8, ...
                     'MarkerEdgeColor',dk, ...
                     'MarkerFaceColor',dk);
        end
    end
end

% now interpolate & plot FOUR lines via PCHIP

nQ = 100;  % points for smooth curve

% 1) Red cluster around 180° (Y in [135,225])
r1 = redY>170 & redY<250;
if sum(r1)>1
  xq = linspace(min(redX(r1)),max(redX(r1)),nQ);
  yq = interp1(redX(r1),redY(r1),xq,'pchip');
  semilogy(xq,yq,'r-','LineWidth',1.2);
end

% 2) Red cluster around 360° (the rest)
r2 = ~r1;
if sum(r2)>1
  xq = linspace(min(redX(r2)),max(redX(r2)),nQ);
  yq = interp1(redX(r2),redY(r2),xq,'pchip');
  semilogy(xq,yq,'r--','LineWidth',1.2);
end

% 3) Green cluster around 90°  (Y<180)
g1 = greenY<180;
if sum(g1)>1
  xq = linspace(min(greenX(g1)),max(greenX(g1)),nQ);
  yq = interp1(greenX(g1),greenY(g1),xq,'pchip');
  semilogy(xq,yq,'g-','LineWidth',1.2);
end

% 4) Green cluster around 270° (Y>=180)
g2 = ~g1;
if sum(g2)>1
  xq = linspace(min(greenX(g2)),max(greenX(g2)),nQ);
  yq = interp1(greenX(g2),greenY(g2),xq,'pchip');
  semilogy(xq,yq,'g--','LineWidth',1.2);
end

xlabel('Date');
ylabel('Deg Error');
hTitle = title('Enceladus Lat/Lon Error & True Anomaly vs Time');

% Increase title (and legend) font sizes
ax     = gca;
baseFS = ax.FontSize;
% Commented out first‐plot legend:
% hLeg   = legend;
% hLeg.FontSize   = 3 * baseFS;
hTitle.FontSize = 3 * baseFS;

% Custom X‐tick labels every 1.37 days (format: day/month)
tickInterval = 1.37;   % days
tickTimes     = Ts(1):tickInterval:Ts(end);
numTicks      = numel(tickTimes);
labels        = strings(1, numTicks);

for i = 1:numTicks
    dv        = mjd20002date(tickTimes(i));    % [Y M D H M S]
    labels(i) = sprintf('%02d/%02d/%02d', dv(3), dv(2), dv(1));  % day/month
end

ax.XTick             = tickTimes;
ax.XTickLabel        = labels;
ax.XTickLabelRotation = 45;
grid on
hold off;

% --- make axis labels & ticks 3× larger for first figure ---
ax = gca;
% baseFS is still the original gca.FontSize before scaling
ax.FontSize       = 3 * baseFS;        % tick labels
ax.XLabel.FontSize = 3 * baseFS;       % x-axis label
ax.YLabel.FontSize = 3 * baseFS;       % y-axis label

%% with true anomaly
figure;
semilogy(Ts, lat_err_series, '-o', ...
         Ts, lon_err_filtered, '-s', ...
         Ts, nu_series, '-^', 'LineWidth', 1.2);
hold on;

% --- add horizontal blue true‐anomaly lines ---
ax     = gca;
baseFS = ax.FontSize;
yline( 90,  'b-', 'LineWidth', 2, 'Label', '90°',  ...
       'FontSize', 2*baseFS, ...
       'LabelHorizontalAlignment','right', ...
       'LabelVerticalAlignment','bottom');
yline(180,  'b-', 'LineWidth', 2, 'Label', '180°', ...
       'FontSize', 2*baseFS, ...
       'LabelHorizontalAlignment','right', ...
       'LabelVerticalAlignment','bottom');
yline(270,  'b-', 'LineWidth', 2, 'Label', '270°', ...
       'FontSize', 2*baseFS, ...
       'LabelHorizontalAlignment','right', ...
       'LabelVerticalAlignment','bottom');
yline(360,  'b-', 'LineWidth', 2, 'Label', '360°', ...
       'FontSize', 2*baseFS, ...
       'LabelHorizontalAlignment','right', ...
       'LabelVerticalAlignment','top');
% ----------------------------------------------

yl = ylim;   % get current y‐limits

% red vertical lines at local minima
for xi = redX
    plot([xi xi], yl, 'r-', 'LineWidth', 1);
end
% green vertical lines at local maxima
for xi = greenX
    plot([xi xi], yl, 'g-', 'LineWidth', 1);
end

% re‐plot the triangles
semilogy(redX,   redY,   'rv', 'MarkerSize', 8, 'MarkerFaceColor','r');
semilogy(greenX, greenY, '^',  'MarkerSize', 8, 'MarkerFaceColor',dk);

xlabel('Date');
ylabel('Deg Errror');
title('Enceladus Lat/Lon Error & True Anomaly (no interp)');
grid on;
hold off;

% --- make axis labels & ticks 3× larger for second figure ---
ax = gca;
baseFS = ax.FontSize;             % capture current font size
ax.FontSize        = 3 * baseFS;  % tick labels
ax.XLabel.FontSize = 3 * baseFS;  % x-axis label
ax.YLabel.FontSize = 3 * baseFS;  % y-axis label