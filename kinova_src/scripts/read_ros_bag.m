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
%     data = readmatrix('ARMOUR_Added_Mass_Reduced_data_v1.csv'); % Reduced Data (Just the Correct Trial)
%     data = readmatrix('ARMOUR_Added_Mass_02_15_2023_v1.csv'); % Raw Data
%     data = readmatrix('ARMOUR_Added_Mass_02_15_2023_v1_Reduced.csv'); % Raw Data
%     data = readmatrix('HardwareDurationValidation03082023.csv');
%     data = readmatrix('ARMOUR_BrakeDebugging_03112023.csv');
%     data = readmatrix('ARMOUR_BrakeDebugging_03112023_12_26_13.csv');
    data = readmatrix('ForceRosBag.csv');

%     processFunction(data)
    compare_desired_v2(data)

end

%% Processing

function output = processFunction(data)

    fig = 0;

    ult_bound_pos = 0.0049; % radians
    ult_bound_vel = 0.09818; % radians / second

    %% Separating Data

    % index offset (for use with raw data since has more columns)
    idx_offset = 3;
    spacing = 1;

    % time
    time = data(:,1);
    % joint position measured?
    joint_pos = data(1:spacing:end,24+idx_offset:30+idx_offset);
    % joint velocity measured?
    joint_vel = data(:,38+idx_offset:44+idx_offset);
    % joint acceleration
    accel = data(:,67+idx_offset:73+idx_offset);
    % joint torque measured
    torque_meas = data(:,2+idx_offset:8+idx_offset);
    % joint torque commanded
    torque_comm = data(:,31+idx_offset:37+idx_offset);

    % joint position error
    joint_pos_err = data(1:spacing:end,9+idx_offset:15+idx_offset);
    % joint velocity error
    joint_vel_err = data(1:spacing:end,16+idx_offset:22+idx_offset);

    %% Processing for NaN values
    for i = 1
        non_nan_idx(:,i) = find(~isnan(joint_pos_err(:,i)));
    end

    %% Plotting

    fig = fig + 1;
    figure(fig)
    hold on
    % plotting joint position error
    plot(time(non_nan_idx), joint_pos_err(non_nan_idx,:))
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
    plot(time(non_nan_idx), joint_vel_err(non_nan_idx,:))
    % plotting velocity error ultimate bound
    plot([time(1) time(end)],[ult_bound_vel ult_bound_vel], '--r')
    plot([time(1) time(end)],[-ult_bound_vel -ult_bound_vel], '--r')
    title('Joint Velocity Error')
    xlabel('Time (s)')
    ylabel('Joint Velocity Error (rad/s)')
    legend('qd1','qd2','qd3','qd4','qd5','qd6','qd7','Velocity Ultimate Bound')
    ylim([-0.15 0.15])

end

function output = compare_desired(data)

    %% Loading Data

    time = data(:,68);
    k_opt = data(:,79:85);
    q_des = data(:,47:53);
    qd_des = data(:,54:60);
    qdd_des = data(:,61:67);

    nan_end = 1;
    time = time(nan_end:end);
    time = time-time(1);
    k_opt = k_opt(nan_end:end,:);
    q_des = q_des(nan_end:end,:);
    qd_des = qd_des(nan_end:end,:);
    qdd_des = qdd_des(nan_end:end,:);

    %% First Iteration
    
    idx_first_iter = 1;
    while norm(k_opt(idx_first_iter,:)) - norm(k_opt(1,:)) == 0
        idx_first_iter = idx_first_iter+1;
    end
    t_first_iter = time(idx_first_iter);

    %% Calculating Desired Trajectory from k_opt

    k_range = [pi/48, pi/48, pi/48, pi/48, pi/48, pi/48, pi/60]'; % needs to match hardware
%     k_range = pi/72*ones(7,1);
%     for j = 1:length(time)
    
        % desired trajectory constraints
        q1 = q_des(idx_first_iter,:)' + k_opt(idx_first_iter,:)' .* k_range; % final position is k*k_range away from initial
        qd1 = zeros(7,1); % final velocity is zero
        qdd1 = zeros(7,1); % final acceleration is zero
    
        duration = 4;
        tid = 100;
        tspan = linspace(0, duration, tid + 1);
        
        beta = match_deg5_bernstein_coefficients({q_des(1,:)', qd_des(1,:)', qdd_des(1,:)', q1, qd1, qdd1}, duration);
    
        for i = 1:length(tspan)
            [q, qd, qdd] = get_desired_traj(beta, tspan(i), duration);
            q_des_matlab(:,i) = q;
            qd_des_matlab(:,i) = qd;
            qdd_des_matlab(:,i) = qdd;
        end
        q_des_matlab = q_des_matlab';
        qd_des_matlab = qd_des_matlab';
        qdd_des_matlab = qdd_des_matlab';

%     end

%     figure(1)
%     title('desired trajectory')
%     subplot(3,1,1)
%     plot(tspan,q_des_matlab)
%     subplot(3,1,2)
%     plot(tspan,qd_des_matlab)
%     subplot(3,1,3)
%     plot(tspan,qdd_des_matlab)

    figure(2)
    hold on
    title('actual trajectory')
    subplot(3,1,1)
    plot(time,q_des)
    subplot(3,1,2)
    plot(time,qd_des)
    subplot(3,1,3)
    plot(time,qdd_des)

%     figure(3)
%     hold on
%     title('Position Comparison First Iteration')
%     plot(tspan,qd_des_matlab)    
%     plot(time(1:idx_first_iter),q_des(1:idx_first_iter,:),'o')
% 
%     figure(4)
%     hold on
%     title('Velocity Comparison First Iteration')
%     plot(tspan,qd_des_matlab)
%     plot(time(1:idx_first_iter),qd_des(1:idx_first_iter,:),'o')
%     % need to fix indices here

    test;


end

function output = compare_desired_v2(data)

    %% Loading Data

    % Debug Topic
    time = data(:,68); % correct index
    debug_q_des = data(:,47:53); % correct indices
    debug_qd_des = data(:,54:60); % correct indices
    debug_qdd_des = data(:,61:67); % correct indices

    % Debug Joint States
    state_pos = data(:,26:32);
    state_vel = data(:,40:46);

    % Debug From Trajectory Topic
    debug_traj_k_opt = data(:,79:85);
    debug_traj_pos = data(:,86:92);
    debug_traj_vel = data(:,96:102);
    debug_traj_accel = data(:,69:75);

    % Debug Joint States Data
    state_pos = data(:,26:32);
    state_vel = data(:,40:46);
    
    % Trajectory Topic
    traj_k_opt = data(:,113:119);
    traj_pos = data(:,120:126);
    traj_vel = data(:,130:136);
    traj_accel = data(:,103:109);

    %% Processing Data

    % Debug Topic Data
    % Offset for NaN Rows
    nan_end = 1;
    % Adjusting Time to Start at Zero
    time = time(nan_end:end);
    time = time-time(1);
    % Adjusting Data to Get Rid of NaN Rows at Beginning
    debug_traj_k_opt = debug_traj_k_opt(nan_end:end,:);
    debug_traj_pos = debug_traj_pos(nan_end:end,:);
    debug_traj_vel = debug_traj_vel(nan_end:end,:);
    debug_traj_accel = debug_traj_accel(nan_end:end,:);
    debug_q_des = debug_q_des(nan_end:end,:);
    debug_qd_des = debug_qd_des(nan_end:end,:);
    debug_qdd_des = debug_qdd_des(nan_end:end,:);
    % Find Other NaN Row Indices
    debug_non_nan_indices = find(~isnan(time));

    % Trajectory Topic Data
    traj_non_nan_idx = find(~isnan(traj_k_opt(:,1)));


    %% First Iteration
    
    for idx_first_iter = 2:length(time)
        if ~(norm(debug_traj_k_opt(idx_first_iter,:)) - norm(debug_traj_k_opt(1,:)) == 0) & ~isnan(debug_traj_k_opt(idx_first_iter,:))
            break
        else
            idx_first_iter = idx_first_iter+1;
        end
    end
    idx_first_iter = idx_first_iter-1;
    t_first_iter = time(idx_first_iter);

    %% Finding Jumps in Desired States

%     j = 1;
%     for i = 2:length(debug_qd_des)
%         debug_q_des_jump(i-1) = debug_q_des(i) - debug_q_des(i-1);
%         debug_qd_des_jump(i-1) = debug_qd_des(i) - debug_qd_des(i-1);
%         debug_qdd_des_jump(i-1) = debug_qdd_des(i) - debug_qdd_des(i-1);
%         if abs(debug_qd_des_jump(i-1)) > 1e-2
%             debug_q_des_jump(i-1)
%             debug_qd_des_jump(i-1)
%             debug_qdd_des_jump(i-1)
%             jump_idx(j) = i-1;
%             j = j + 1;
%         end
%     end
% 
%     brake_indices = jump_idx;

    %% Calculating Desired Trajectory from k_opt

    k_range = [pi/48, pi/48, pi/48, pi/48, pi/48, pi/48, pi/60]'; % needs to match hardware
%     k_range = pi/72*ones(7,1);
    
    % desired trajectory constraints
    q1 = debug_q_des(idx_first_iter-1,:)' + debug_traj_k_opt(idx_first_iter-1,:)' .* k_range; % final position is k*k_range away from initial
    qd1 = zeros(7,1); % final velocity is zero
    qdd1 = zeros(7,1); % final acceleration is zero

    duration = 1.5;
    tid = 100;
    tspan = linspace(0, duration, tid + 1);
    
    beta = match_deg5_bernstein_coefficients({debug_q_des(1,:)', debug_qd_des(1,:)', debug_qdd_des(1,:)', q1, qd1, qdd1}, duration);

    for i = 1:length(tspan)
        [q, qd, qdd] = get_desired_traj(beta, tspan(i), duration);
        q_des_matlab(:,i) = q;
        qd_des_matlab(:,i) = qd;
        qdd_des_matlab(:,i) = qdd;
    end
    q_des_matlab = q_des_matlab';
    qd_des_matlab = qd_des_matlab';
    qdd_des_matlab = qdd_des_matlab';

    %% Plotting 

    % Desired Variables
%     figure(1)
%     hold on
%     title('debug desired')
%     % Desired Position
%     subplot(3,1,1)
%     hold on
%     plot(time,debug_q_des)
% %     plot([time(brake_indices) time(brake_indices)],[-5 5],'--r')
%     % Desired Velocity
%     subplot(3,1,2)
%     hold on
%     plot(time,debug_qd_des)
% %     plot([time(brake_indices) time(brake_indices)],[-0.15 0.15],'--r')
%     % Desired Acceleration
%     subplot(3,1,3)
%     hold on
%     plot(time,debug_qdd_des)
% %     plot([time(brake_indices) time(brake_indices)],[-0.5 0.5],'--r')

    % Plot Measured Joint States
%     figure(2)
%     title('Joint States')
%     subplot(2,1,1)
%     plot(time,state_pos)
%     subplot(2,1,2)
%     plot(time,state_vel)

    % Plotting k_opt
%     figure(3)
%     hold on
%     title('k_{opt}')
%     plot(time,debug_traj_k_opt)
%     plot(time(traj_non_nan_idx-1)+1.5,traj_k_opt(traj_non_nan_idx,:),'--o')
    % Looks like all the k_opt are published before they are needed and the
    % ones in debug match up in value.

    % Plot Trajectory Message Velocity vs Desired Velocity from Debug
%     figure(4)
%     hold on
%     title('Desired Velocity Comparison')
%     plot(time,debug_qd_des)
%     plot(time(traj_non_nan_idx-1)+1.5,traj_vel(traj_non_nan_idx,:),'--o')

%     figure(3)
%     title('k_{opt}')
%     hold on
%     plot(time,debug_k_opt)
%     plot(time(non_nan_idx(1:end)-1),debug_traj_k_opt(non_nan_idx(1:end),:),'o')
% 
    figure(6)
    title('Desired Velocity Comparison')
    hold on
    plot(time(1:idx_first_iter,:),debug_qd_des(1:idx_first_iter,:))
    plot(tspan, qd_des_matlab,'o')

    test;


end


%% helper functions
function [q, qd, qdd] = get_desired_traj(beta, t, duration)

    if nargin < 3
        duration = 1;
    end

    [B, dB, ddB] = Bezier_kernel_deg5(t/duration); %t/dur

    q = zeros(7,length(t));
    qd = zeros(7,length(t));
    qdd = zeros(7,length(t));
    
    for j = 1:6
        q = q + beta{j} .* B(:,j)';
        qd = qd + beta{j} .* dB(:,j)';
        qdd = qdd + beta{j} .* ddB(:,j)';
    end

    qd = qd / duration; % commenting out these durations matches the c++ and matlab code
    qdd = qdd / duration / duration;
end