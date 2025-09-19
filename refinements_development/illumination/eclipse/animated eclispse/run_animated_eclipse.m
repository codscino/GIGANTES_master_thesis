
clc; close all; clear all;

% script to run animated scatter plot of Saturn eclipsing
% Enceladus eclipse at equinox

start_date = [2025, 5, 10, 21, 46, 30];
end_date = [2025, 5, 10, 21, 47, 40];

T_start = date2mjd2000(start_date);
T_end = date2mjd2000(end_date);

time_step = 1; 
grid_angle_step = 0.5; % put 2 for faster plotting

parallel = true; % use parallel computing

animated_eclipse(T_start, T_end, time_step, grid_angle_step, parallel);