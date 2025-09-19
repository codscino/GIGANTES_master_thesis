clc
clear

%% decide time
date1 = [2025, 5, 5, 12, 0, 0];
date2 = [2125, 5, 5, 12, 0, 0];
jd_date1 = date2mjd2000(date1);
jd_date2 = date2mjd2000(date2);

%% load kernels
id = 699; % Saturn
% Time check
spiceCheck([jd_date1,jd_date2], id)

% compute position
[r_mine,v_mine,kep] = EphSS_car_spice(id,jd_date1,true);

%% mice direct call
cspice_kclear;    
SPKfile           = 'sat441.bsp';
leapsecondsKernel = 'naif0012.tls';
kernels           = {SPKfile, leapsecondsKernel};
spkPath           = which(SPKfile);
lskPath           = which(leapsecondsKernel);
cspice_furnsh({ lskPath, spkPath });

spiceParam.frame    = 'J2000';
spiceParam.abcorr   = 'NONE';
spiceParam.observer = '0';  % SSB

spice_time = jd_date1*86400;

[state, ~] = cspice_spkezr( ...
    num2str(id), ...
    spice_time, ...
    spiceParam.frame, ...
    spiceParam.abcorr, ...
    spiceParam.observer );
state = state';
r_mice = state(1:3);
v_mice = state(4:6);

%% error
r_error = abs(r_mine-r_mice)
v_error = abs(v_mine-v_mice)

%% for loop for speed
% number of daily steps (inclusive)
% nDays = round(jd_date2 - jd_date1) + 1;
nDays = 3652500;

% pre-allocate: 3 rows (x,y,z) × nDays columns
r = zeros(3, nDays);
v = zeros(3, nDays);

% loop with an explicit index
tic
for i = 1:nDays
    jd = jd_date1 + (i-1)/100;
    [r(:,i), v(:,i)] = EphSS_car_spice(id, jd, true);
end
toc %0.35 seconds on Claudio Mac

