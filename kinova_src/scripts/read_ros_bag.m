%% AMROUR ROS Bag Reader
% Zachary Brei
% 02/13/2023

clear all; 
close all; clc;

%% Loading the data

try

    test1 = data(1);

    processFunction(data);

catch

    warning('variable data not loaded yet, loading now')

%     data = readmatrix('ARMOUR_Added_Mass_v1.csv'); % Raw Data
    data = readmatrix('ARMOUR_Added_Mass_Reduced_data_v1.csv'); % Reduced Data (Just the Correct Trial)

    processFunction(data)

end

%% Processing

function output = processFunction(data)

    fig = 0;

    ult_bound_pos = 0.0049; % radians
    ult_bound_vel = 0.09818; % radians / second

    %% Separating Data

    % index offset (for use with raw data since has more columns)
    idx_offset = 3;
    
    % time
    time = data(:,1);
    % joint position measured?
    joint_pos = data(:,24+idx_offset:30+idx_offset);
    % joint velocity measured?
    joint_vel = data(:,38+idx_offset:44+idx_offset);
    % joint acceleration
    accel = data(:,67+idx_offset:73+idx_offset);
    % joint torque measured
    torque_meas = data(:,2+idx_offset:8+idx_offset);
    % joint torque commanded
    torque_comm = data(:,31+idx_offset:37+idx_offset);

    % joint position error
    joint_pos_err = data(:,9+idx_offset:15+idx_offset);
    % joint velocity error
    joint_vel_err = data(:,16+idx_offset:22+idx_offset);

    %% Plotting

    fig = fig + 1;
    figure(fig)
    hold on
    % plotting joint position error
    plot(time, joint_pos_err)
    % plotting position error ultimate bound
    plot([time(1) time(end)],[ult_bound_pos ult_bound_pos], '--r')
    plot([time(1) time(end)],[-ult_bound_pos -ult_bound_pos], '--r')
    title('Joint Position Error')
    xlabel('Time (s)')
    ylabel('Joint Position Error (rad)')
    legend('q1','q2','q3','q4','q5','q6','q7','Position Ultimate Bound')
    ylim([-0.0075 0.0075])

    fig = fig + 1;
    figure(fig)
    hold on
    % plotting joint velocity error
    plot(time, joint_vel_err)
    % plotting velocity error ultimate bound
    plot([time(1) time(end)],[ult_bound_vel ult_bound_vel], '--r')
    plot([time(1) time(end)],[-ult_bound_vel -ult_bound_vel], '--r')
    title('Joint Velocity Error')
    xlabel('Time (s)')
    ylabel('Joint Velocity Error (rad/s)')
    legend('qd1','qd2','qd3','qd4','qd5','qd6','qd7','Velocity Ultimate Bound')
    ylim([-0.15 0.15])

end