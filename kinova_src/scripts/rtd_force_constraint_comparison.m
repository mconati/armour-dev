%% RTD-Force Matlab and Realtime Constraint Comparison
% Zachary Brei
% 01/13/2023

clear all; close all; clc;

%% Robot Parameters
agent_urdf = 'Kinova_Grasp_URDF.urdf';

add_uncertainty_to = 'all'; % choose 'all', 'link', or 'none'
links_with_uncertainty = {}; % if add_uncertainty_to = 'link', specify links here.
uncertain_mass_range = [0.97, 1.03];

robot = importrobot(agent_urdf);
robot.DataFormat = 'col';
robot.Gravity = [0 0 -9.81];
params = load_robot_params(robot, ...
                           'add_uncertainty_to', add_uncertainty_to, ...
                           'links_with_uncertainty', links_with_uncertainty,...
                           'uncertain_mass_range', uncertain_mass_range);

for i=1:7
    % lower limit
    joint_position_limits(1,i) = robot.Bodies{1, i}.Joint.PositionLimits(1);
    % upper limit
    joint_position_limits(2,i) = robot.Bodies{1, i}.Joint.PositionLimits(2);
end

joint_speed_limits = [-1.3963, -1.3963, -1.3963, -1.3963, -1.2218, -1.2218, -1.2218;
                       1.3963,  1.3963,  1.3963,  1.3963,  1.2218,  1.2218,  1.2218]; % matlab doesn't import these from urdf so hard code into class
joint_input_limits = [-56.7, -56.7, -56.7, -56.7, -29.4, -29.4, -29.4;
                       56.7,  56.7,  56.7,  56.7,  29.4,  29.4,  29.4]; % matlab doesn't import these from urdf so hard code into class
transmision_inertia = [8.02999999999999936 11.99620246153036440 9.00254278617515169 11.58064393167063599 8.46650409179141228 8.85370693737424297 8.85873036646853151]; % matlab doesn't import these from urdf so hard code into class
M_min_eigenvalue = 5.095620491878957; % matlab doesn't import these from urdf so hard code into class

max_travel = 0.4;

% obstacles = {};
obstacles{1} = box_obstacle_zonotope('center', [3; 0; 0.6],...
                                     'side_lengths', [0.1; 0.1; 0.1]) ;

goal_type = 'configuration'; % pick 'end_effector_location' or 'configuration'
goal_radius = pi/30;
dimension = 3 ;
verbosity = 10;

% P = uarmtd_planner('verbose', verbosity, ...
%                    'first_iter_pause_flag', true, ...
%                    'use_q_plan_for_cost', true, ...
%                    'input_constraints_flag', true, ...
%                    'use_robust_input', use_robust_input, ...
%                    'traj_type', 'bernstein', ...
%                    'use_cuda', false) ;

% W = kinova_world_static('create_random_obstacles_flag', false, 'goal_radius', goal_radius, 'N_obstacles',length(obstacles),'dimension',dimension,'workspace_goal_check', 0,...
%                         'verbose',verbosity, 'start', start, 'goal', goal, 'obstacles', obstacles, 'goal_type', goal_type) ;


%% Create list of random conditions to test

num_conditions = 1;

for i = 1:num_conditions
    % r = a + (b-a).*rand(N,1)
    for j = 1:7
        
        if joint_position_limits(2,j) > 100
            % random start positions
            q_0(i,j) = rand();
            % random goal positions
            q_goal(i,j) = rand();
        else
            % random start positions
            q_0(i,j) = joint_position_limits(1,j) + (joint_position_limits(2,j)-joint_position_limits(1,j))*rand();
            % random goal positions
            q_goal(i,j) = joint_position_limits(1,j) + (joint_position_limits(2,j)-joint_position_limits(1,j))*rand();
        end
        
        % generate waypoint for desired position
        dir = q_goal(i,j)-q_0(i,j);
        dir = dir./norm(dir);
        q_des(i,j) = q_0(i,j)+max_travel*dir;
        % random start velocity
        q_dot_0(i,j) = joint_speed_limits(1,j) + (joint_speed_limits(2,j)-joint_speed_limits(1,j))*rand(1,1);
        % random start acceleration (what interval?)
        q_ddot_0(i,j) = -1 + (1+1)*rand(1,1);
    end
end

%% C++ Method

%% Calling Realtime Planner

% !!!!!!
% needs to match c++ code
% !!!!!!
% P.jrs_info.n_t = 128;
% P.jrs_info.n_q = 7;
% P.jrs_info.n_k = 7;
% P.jrs_info.c_k_bernstein = zeros(7,1);

% !!!!!!
% Make sure this is consistent with the k_range in
% cuda-dev/PZsparse-Bernstein/Trajectory.h 
% !!!!!!
% P.jrs_info.g_k_bernstein = [pi/72; pi/72; pi/72; pi/72; pi/72; pi/72; pi/72];

for idx_real = 1:num_conditions

    % organize input to cuda program
    fprintf('Calling CUDA & C++ Program!')
    cuda_input_file = fopen('/home/roahmlab/Documents/armour-dev/kinova_src/kinova_simulator_interfaces/kinova_planner_realtime/buffer/armour.in', 'w');
    
    for ind = 1:length(q_0)
        fprintf(cuda_input_file, '%.10f ', q_0(idx_real,ind));
    end
    fprintf(cuda_input_file, '\n');
    for ind = 1:length(q_dot_0)
        fprintf(cuda_input_file, '%.10f ', q_dot_0(idx_real,ind));
    end
    fprintf(cuda_input_file, '\n');
    for ind = 1:length(q_ddot_0)
        fprintf(cuda_input_file, '%.10f ', q_ddot_0(idx_real,ind));
    end
    fprintf(cuda_input_file, '\n');
    for ind = 1:length(q_des)
        fprintf(cuda_input_file, '%.10f ', q_des(idx_real,ind));
    end
    fprintf(cuda_input_file, '\n');
    fprintf(cuda_input_file, '%d\n', max(length(obstacles), 0));
    for obs_ind = 1:length(obstacles)
        temp = reshape(obstacles{obs_ind}.Z, [1,size(obstacles{obs_ind}.Z,1) * size(obstacles{obs_ind}.Z,2)]); %world_info.
        for ind = 1:length(temp)
            fprintf(cuda_input_file, '%.10f ', temp(ind));
        end
        fprintf(cuda_input_file, '\n');
    end
    
    fclose(cuda_input_file);
    
    % call cuda program in terminal
    % you have to be in the proper path!
    %                     terminal_output = system('./../kinova_simulator_interfaces/kinova_planner_realtime/armour_main'); % armour path
    terminal_output = system('/home/roahmlab/Documents/armour-dev/kinova_src/kinova_simulator_interfaces/kinova_planner_realtime/rtd_force_main'); % rtd-force path
    
    % To do, figure out how to read and store output for comparison
    
    if terminal_output == 0
        data = readmatrix('armour.out', 'FileType', 'text');
        k_opt = data(1:end-1);
        planning_time(i) = data(end) / 1000.0; % original data is milliseconds
    
        if length(k_opt) == 1
            fprintf('Unable to find new trajectory!')
            k_opt = nan;
        else
            fprintf('New trajectory found!');
            for i = 1:length(k_opt)
                fprintf('%7.6f ', k_opt);
            end
            fprintf('\n');
        end
    else
        error('CUDA program error! Check the executable path in armour-dev/kinova_src/kinova_simulator_interfaces/uarmtd_planner');
    end
    
    if terminal_output == 0
        % read FRS information if needed
        joint_frs_center(:,i) = readmatrix('armour_joint_position_center.out', 'FileType', 'text');
        joint_frs_radius(:,i) = readmatrix('armour_joint_position_radius.out', 'FileType', 'text');
        control_input_radius(:,i) = readmatrix('armour_control_input_radius.out', 'FileType', 'text');
        constraints_value(:,i) = readmatrix('armour_constraints.out', 'FileType', 'text');
    else
        k_opt = nan;
    end
    k_opt_storage(:,i) = k_opt;
end


%% Matlab Method


% generate constraints

for i=1:num_conditions

    matlab_constraints{i} = generate_constraints(P, q_0, q_dot_0, q_ddot_0, O, q_des);

end

