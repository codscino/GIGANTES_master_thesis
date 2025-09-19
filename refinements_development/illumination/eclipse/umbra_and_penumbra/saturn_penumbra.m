clear
clc
r_sun = getAstroConstants('Sun', 'radius');

%% Earth - Moon
% geosynchronous
r_earth = getAstroConstants('Earth', 'radius');
au = getAstroConstants('AU');
d_eg = 35786; %earth-geosynchronous distance [km]   d_eg = 42163 - r_earth
[~, h_pen_geo] = umbr_pen(r_sun, r_earth, au, d_eg); 
vallado = 13098;
error = 13098 - h_pen_geo*2;

% eclipse
d_em = 384400; % earth-moon distance [km]
[h_umb1, h_pen1] = umbr_pen(r_sun, r_earth, au, d_em);
ratio1 = h_umb1/h_pen1;

%% Saturn - Enceladus
% constants
r_saturn = getAstroConstants('Saturn', 'radius');
a_saturn = getAstroConstants('Saturn', 'Sma');
d_se = 238200; 

% eclipse
[h_umb2, h_pen2] = umbr_pen(r_sun, r_saturn, a_saturn, d_se);
ratio2 = h_umb2/h_pen2;


