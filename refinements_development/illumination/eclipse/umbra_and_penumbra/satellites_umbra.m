% As script to calculat the umbra and penumbra ratios of the eclipse casted
% by Saturn and it moons on Enceladus.
% The umb_norm is the h_umb normalize with the % Enceladus radius: 252.100 km
%
% Saturn    : h_umb = 58125.170 km, umb_norm = 230.564, umb/pen ratio =  0.996
% Mimas     : h_umb = 107.486 km, umb_norm =  0.426, umb/pen ratio =  0.372
% Tethys    : h_umb = 387.050 km, umb_norm =  1.535, umb/pen ratio =  0.573
% Dione     : h_umb = 376.851 km, umb_norm =  1.495, umb/pen ratio =  0.505
% Rhea      : h_umb = 506.066 km, umb_norm =  2.007, umb/pen ratio =  0.495
% Titan     : h_umb = 1978.642 km, umb_norm =  7.849, umb/pen ratio =  0.623
% 
% --- Moon(shadowed by Earth):  umb/pen Ratio =   0.561

clear
clc

r_sun       = getAstroConstants('Sun',       'radius');
r_saturn    = getAstroConstants('Saturn',    'radius');
a_saturn    = getAstroConstants('Saturn',    'Sma');

%% Saturn umbra
[d_sat_enceladus, ~, r_enceladus, ~] = satMoonsConstants(1);
[h_umb_saturn, h_pen_saturn] = umbr_pen(r_sun, r_saturn, a_saturn, d_sat_enceladus);
ratio_sat = h_umb_saturn / r_enceladus;

% Print Enceladus radius
fprintf('Enceladus radius: %.3f km\n', r_enceladus);

% Calculate umb/pen ratio for Saturn
umb_pen_ratio_saturn = h_umb_saturn / h_pen_saturn;

% Now modify Saturn’s shadow stats print line
fprintf('%-10s: h_umb = %7.3f km, umb_norm = %6.3f, umb/pen ratio = %6.3f\n', ...
        'Saturn', h_umb_saturn, ratio_sat, umb_pen_ratio_saturn);

%% for‐loop moons
% d_sm → saturn-moon distance
% r_m  → radius of the moon

% Define the moon names and IDs
moonNames = {'Mimas','Enceladus','Tethys','Dione','Rhea','Titan'};
moonIDs   = 0:5;            % [0=Mimas, 1=Enceladus, …, 5=Titan]

% Remove Enceladus (ID == 1) in one line:
mask      = moonIDs~=1;    
moonNames = moonNames(mask);
moonIDs   = moonIDs(mask);

nMoons    = numel(moonIDs);
h_umb     = zeros(1,nMoons);
h_pen     = zeros(1,nMoons);
h_umb_norm= zeros(1,nMoons);
umb_pen_ratio = zeros(1,nMoons);

for k = 1:nMoons
    id = moonIDs(k);
    [d_sm, ~, r_m, ~] = satMoonsConstants(id);
    [h_umb(k), h_pen(k)] = umbr_pen(r_sun, r_m, a_saturn, abs(d_sm - r_enceladus));
    h_umb_norm(k)   = h_umb(k) / r_enceladus;    % normalized to Enceladus radius
    umb_pen_ratio(k)= h_umb(k) / h_pen(k);       % umb/pen ratio
end

% Display
for k = 1:nMoons
    fprintf('%-10s: h_umb = %7.3f km, umb_norm = %6.3f, umb/pen ratio = %6.3f\n', ...
            moonNames{k}, h_umb(k), h_umb_norm(k), umb_pen_ratio(k));
end

%% Earth - Moon
% eclipse
d_em = 384400; % earth-moon distance [km]
r_earth = getAstroConstants('Earth', 'radius'); % Ensure r_earth is defined if not already from previous sections
au = getAstroConstants('AU'); % Ensure au is defined

[h_umb1, h_pen1] = umbr_pen(r_sun, r_earth, au, d_em);
ratio1 = h_umb1/h_pen1;

fprintf('\n--- Moon(shadowed by Earth):');
fprintf('  umb/pen Ratio = %7.3f\n', ratio1);
