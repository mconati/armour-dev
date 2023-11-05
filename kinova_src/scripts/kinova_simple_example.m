%% description
% This script defines a simple world with user defined obstacles and runs
% the ARMOUR planner on them. It then saves information on how well each the
% planner performed in each trial.
%
% Authors: Bohao Zhang (adapted from Patrick Holmes code)
% Created 25 November 2022

initialize_script_path = matlab.desktop.editor.getActiveFilename;
cd(initialize_script_path(1:end-23));

close all; clear; clc;

%% user parameters
goal_type = 'configuration'; % pick 'end_effector_location' or 'configuration'
goal_radius = pi/30;
dimension = 3 ;
verbosity = 10;

%%% for planner
traj_type = 'bernstein'; % pick 'orig' (ARMTD) or 'bernstein' (ARMOUR)
use_cuda_flag = true;

%%% for agent
agent_urdf = 'kinova_without_gripper.urdf';

add_uncertainty_to = 'all'; % choose 'all', 'link', or 'none'
links_with_uncertainty = {}; % if add_uncertainty_to = 'link', specify links here.
uncertain_mass_range = [0.97, 1.03];

agent_move_mode = 'save_full_trajectories' ; % pick 'direct' or 'integrator' or 'save_full_trajectories'
use_CAD_flag = true; % plot robot with CAD or bounding boxes

%%% for LLC
use_robust_input = true;
LLC_V_max = 5e-5;

%%% for HLP
if_use_RRT = false;
HLP_grow_tree_mode = 'new' ; % pick 'new' or 'keep'
plot_waypoint_flag = true ;
plot_waypoint_arm_flag  = true ;
lookahead_distance = 0.1 ;

% plotting
plot_while_running = false ;

% simulation
max_sim_time = 172800 ; % 48 hours
max_sim_iter = 600 ;
stop_threshold = 4 ; % number of failed iterations before exiting

%%% for world
start = [0; -pi/2; 0; 0; 0; 0; 0]; % start configuration
goal = [pi/4; -pi/2; 0; 0; 0; 0; 0]; % goal configuration


obstacles{1} = box_obstacle_zonotope('center', [0; 0; 0.6],...
                                     'side_lengths', [0.1; 0.1; 0.1]) ;
obstacles{2} = box_obstacle_zonotope('center', [0.3; 0; 0.4],...
                                     'side_lengths', [0.1; 0.8; 0.05]) ;

%% robot params:
robot = importrobot(agent_urdf);
robot.DataFormat = 'col';
robot.Gravity = [0 0 -9.81];
params = load_robot_params(robot, ...
                           'add_uncertainty_to', add_uncertainty_to, ...
                           'links_with_uncertainty', links_with_uncertainty,...
                           'uncertain_mass_range', uncertain_mass_range);
joint_speed_limits = [-1.3963, -1.3963, -1.3963, -1.3963, -1.2218, -1.2218, -1.2218;
                       1.3963,  1.3963,  1.3963,  1.3963,  1.2218,  1.2218,  1.2218]; % matlab doesn't import these from urdf so hard code into class
joint_input_limits = [-56.7, -56.7, -56.7, -56.7, -29.4, -29.4, -29.4;
                       56.7,  56.7,  56.7,  56.7,  29.4,  29.4,  29.4]; % matlab doesn't import these from urdf so hard code into class
transmision_inertia = [8.02999999999999936 11.99620246153036440 9.00254278617515169 11.58064393167063599 8.46650409179141228 8.85370693737424297 8.85873036646853151]; % matlab doesn't import these from urdf so hard code into class
M_min_eigenvalue = 5.095620491878957; % matlab doesn't import these from urdf so hard code into class

%% automated from here
if plot_while_running
    figure(1); clf; view(3); grid on;
end

% run loop
tic;
W = kinova_world_static('create_random_obstacles_flag', false, 'goal_radius', goal_radius, 'N_obstacles',length(obstacles),'dimension',dimension,'workspace_goal_check', 0,...
                        'verbose',verbosity, 'start', start, 'goal', goal, 'obstacles', obstacles, 'goal_type', goal_type) ;

% create arm agent
A = uarmtd_agent(robot, params,...
                 'verbose', verbosity,...
                 'animation_set_axes_flag', 0,... 
                 'animation_set_view_flag', 0,...
                 'move_mode', agent_move_mode,...
                 'use_CAD_flag', use_CAD_flag,...
                 'joint_speed_limits', joint_speed_limits, ...
                 'joint_input_limits', joint_input_limits, ...
                 'add_measurement_noise_', false, ...
                 'measurement_noise_size_', 0,...
                 'M_min_eigenvalue', M_min_eigenvalue, ...
                 'transmision_inertia', transmision_inertia);

% LLC
if use_robust_input
    A.LLC = uarmtd_robust_CBF_LLC('verbose', verbosity, ...
                                  'use_true_params_for_robust', false, ...
                                  'V_max', LLC_V_max, ...
                                  'if_use_mex_controller', true);
else
    A.LLC = uarmtd_nominal_passivity_LLC('verbose', verbosity);
end

A.LLC.setup(A);

P = uarmtd_planner('verbose', verbosity, ...
                   'first_iter_pause_flag', false, ...
                   'use_q_plan_for_cost', true, ...
                   'input_constraints_flag', true, ...
                   'use_robust_input', use_robust_input, ...
                   'traj_type', traj_type, ...
                   'use_cuda', use_cuda_flag, ...
                   'waypoint_mode', 'goal', ...
                   'save_constraints', true) ;

if if_use_RRT
    P.HLP = arm_end_effector_RRT_star_HLP('plot_waypoint_flag',plot_waypoint_flag,...
                                          'plot_waypoint_arm_flag',plot_waypoint_arm_flag,...
                                          'grow_tree_mode',HLP_grow_tree_mode,...
                                          'buffer',0.1) ;
end

% set up world using arm
I = A.get_agent_info ;
W.setup(I) ;
W.bounds = [-1 1 -1 1 0 2];

% place arm at starting configuration
A.state(A.joint_state_indices) = W.start ;

% create simulator
S = simulator_armtd(A,W,P, ...
                    'verbose', verbosity, ...
                    'stop_threshold', stop_threshold, ...
                    'plot_while_running', plot_while_running,...
                    'allow_replan_errors',true,...
                    'max_sim_time',max_sim_time,...
                    'max_sim_iterations',max_sim_iter,...
                    'stop_sim_when_ultimate_bound_exceeded', use_robust_input) ; 

% %% plotting
if plot_while_running
    figure(1) ; clf ; axis equal ; xlim([-1 1]); ylim([-1 1]); zlim([0 2]); grid on; hold on ;

    if dimension == 3
        view(3);
    end
    
    plot(A);
    plot(W);
end

% run simulation
summary = S.run();
%% Plotting States

figure(101)
subplot(3,1,1)
title('Joint Positions')
plot(A.time,A.state(A.joint_state_indices,:))
subplot(3,1,2)
title('Joint Velocities')
plot(A.time,A.state(A.joint_speed_indices,:))
subplot(3,1,3)
title('Joint Torques')
plot(A.time,A.input)

%% 

% Calculating the Acceleration

joint_angles = A.full_state(A.joint_state_indices,:);
joint_angular_velocity = A.full_state(A.joint_speed_indices,:);

qdd_post = zeros(7,length(A.full_time));
% calculating the acceleration in post to compare with what is stored
for i = 2:length(A.full_time)
    [M, C, g] = A.calculate_dynamics(joint_angles(:,i), joint_angular_velocity(:,i), A.params.true);

    for j = 1:A.n_inputs
        M(j,j) = M(j,j) + A.transmision_inertia(j);
    end
    % can I call u=A.LLC.get_control_inputs() here with the P.info?
    
    qdd_post(:,i) = M\(A.full_u(:,i)-C*joint_angular_velocity(:,i)-g);
%         qdd_post(:,i) = M\(A.input(:,i)-C*joint_angular_velocity(:,i)-g);

end


%% Plotting Full trajectories
%states
T = A.full_time;
Ys = A.full_state(A.joint_state_indices, :);

%velocities
Yv = A.full_state(A.joint_speed_indices, :);

%accelerations
Ya = qdd_post;

%inputs
%Yu = A.full_u;

figure(102)
subplot(3,1,1)
title('Joint Positions')
plot_whole_trajectories_with_bounds(A, T ,Ys, A.position_bounds(:, 1), false);
subplot(3,1,2)
title('Joint Velocities')
plot_whole_trajectories_with_bounds(A, T ,Yv, A.velocity_bounds(:, 1), false);
subplot(3,1,3)
title('Joint Accelerations')
plot_whole_trajectories_with_bounds(A, T ,Ya, -1, false);

%plot torques for the first trajectory
%plot_torques(A.full_u, A.input_constraints, A.input_radii, A.full_time, 1, 0)

%Plot for all trajectories
%plot_torques(A.full_u, A.input_constraints, A.input_radii, A.full_time, -1, 0)

% %%
% function plot_whole_trajectories(A, T, Y, makefigure)
%     buffer_y = [Y(:, 1)];
%     buffer_yt = [T(:, 1)];
%     buffer_b = [];
%     buffer_bt = [];
%     t = A.time;
%     
%     t_stop = A.t_stop; 
%     
%     % Create a figure and axis
%     if makefigure
%        figure;
%     end
%     hold on;
%     
%     color_idx = 1; % Initialize the color index
%     start_time = 0;
%     s_T = 1;
%     while start_time < max(T)
%         end_time = start_time + t_stop;
%         % Clip the time interval to avoid going beyond the data
%         if end_time > max(T)
%             end_time = max(T);
%         end
%         indices = find(T > start_time & T <= end_time);
%         % Select the appropriate color
%         if color_idx == 1
%             color = 'b'; % Blue
%             disp(indices)
%             buffer_y = [buffer_y, Y(:, indices)];
%             buffer_yt = [buffer_yt, T(:, s_T:(s_T+length(indices))-1)];
%             s_T = s_T + length(indices)-1;
%         else
%             color = 'r'; % Red
%             buffer_b = [buffer_b, Y(:, indices)];
%             buffer_bt = [buffer_bt, T(:, s_T:(s_T+length(indices))-1)];
%             plot(T(:, s_T:(s_T+length(indices)-1)), Y(:, indices), color);
%         end
%         
%         % Find the corresponding indices in the time array
%     
%         
%         % Plot the line segment
%     
%         
%         % Update the color index and start time for the next segment
%         color_idx = 3 - color_idx; % Toggle between 1 and 2
%         start_time = end_time;
%     
%     end
%     plot(buffer_yt, buffer_y, 'b');
%     xlabel('Time (seconds)');
%     ylabel('Y');
%     grid on;
% end
% 
% %% 
%  