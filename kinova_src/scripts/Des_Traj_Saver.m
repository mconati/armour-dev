%% Saving Evenly Sampled Desired Trajectory for Python Rendering

clear all;
close all;
clc;

%% Load Data

load('HardwareVideo3Runs_v2.mat')

%% Extract the Data

time = rosbag.raw_time - rosbag.raw_time(1);
q_des = (rosbag.debug_q_des);

%% Interpolate the Data

framerate = 50;
test = linspace(0,0.875,framerate);
spacing = diff(test);
% time(end) = 
interp_time = 0:spacing(1):0.875*395; % linspace(0,0.875*395,spacing(1));

interp_q_des = interp1(time,q_des,interp_time);

%% Plotting to Check

figure(1)
hold on
plot(time,q_des)
plot(interp_time,interp_q_des,'or')

%% Saving Data

% python expects nx7 array
filename = 'HardwareVideo_05_05_2023_q_des.csv';
writematrix(interp_q_des,filename,'Delimiter',',');