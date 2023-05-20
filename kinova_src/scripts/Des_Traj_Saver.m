%% Saving Evenly Sampled Desired Trajectory for Python Rendering

clear all;
close all;
clc;

%% Load Data

load('HardwareVideo_MultipleTrials_05_17_2023_ROSData.mat')

%% Extract the Data

time = rosbag.raw_time - rosbag.raw_time(1);
q_des = (rosbag.debug_q_des);

%% Interpolate the Data

framerate = 50;
test = linspace(0,0.875,framerate+1);
spacing = diff(test);
% time(end) = 
interp_time = 0:spacing(1):0.875*332;
% Note that the number that multiplies the planning time is calculated by
% dividing the total time of the trials by the planning time i.e.:
% time(end) / planning time

interp_q_des = interp1(time,q_des,interp_time);

%% Plotting to Check

figure(1)
hold on
plot(time,q_des)
plot(interp_time,interp_q_des,'or')

%% Saving Data

% python expects nx7 array
filename = 'HardwareVideo_MultipleTrials_05_17_2023_ROSData_q_des.csv';
% writematrix(interp_q_des,filename,'Delimiter',',');