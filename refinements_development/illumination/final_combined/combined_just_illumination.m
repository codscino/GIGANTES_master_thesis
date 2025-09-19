%% plot_illumination
clc; clear all; close all;

% init
date = [2025, 5, 10, 21, 46, 59];
T = date2mjd2000(date);
kernels = {'sat441.bsp', 'naif0012.tls'};
loadSpiceKernels(kernels);

parallel = true;
step_deg = 0.5;
FF = 27; % terminator+eclipse

plot_illumination(T, FF, step_deg, parallel, true)