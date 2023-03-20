%% description
% This script defines a simple world with user defined obstacles and runs
% the ARMOUR planner on them. It then saves information on how well each the
% planner performed in each trial.
%
% Authors: Bohao Zhang (adapted from Patrick Holmes code)
% Created 25 November 2022

initialize_script_path = matlab.desktop.editor.getActiveFilename;
cd(initialize_script_path(1:end-29));

close all; clear; clc; dbstop if error

%% user parameters
goal_type = 'configuration'; % pick 'end_effector_location' or 'configuration'
goal_radius = pi/30;
dimension = 3 ;
verbosity = 10;

DURATION = 1.0;
k_range = pi/72; % used for graph planner check

u_s = 0.609382421; 
surf_rad =  0.058 / 2;

%%% for planner
traj_type = 'bernstein'; % pick 'orig' (ARMTD) or 'bernstein' (ARMOUR)
use_cuda_flag = false;

%%% for agent
agent_urdf = 'Kinova_Grasp_URDF.urdf';

add_uncertainty_to = 'all'; % choose 'all', 'link', or 'none'
links_with_uncertainty = {}; % if add_uncertainty_to = 'link', specify links here.
uncertain_mass_range = [0.97, 1.03];

agent_move_mode = 'integrator' ; % pick 'direct' or 'integrator'
use_CAD_flag = true; % plot robot with CAD or bounding boxes

%%% for LLC
use_robust_input = false;
LLC_V_max = 1e-2;
alpha_constant = 10;
Kr = 5;

%%% for HLP
if_use_RRT = false;
HLP_grow_tree_mode = 'new' ; % pick 'new' or 'keep'
plot_waypoint_flag = true ;
plot_waypoint_arm_flag  = true ;
lookahead_distance = 0.1 ;

if_use_graph_planner = true;

% plotting
plot_while_running = true ;

% simulation
max_sim_time = 172800 ; % 48 hours
max_sim_iter = 600 ;
stop_threshold = 4 ; % number of failed iterations before exiting

obstacles{1} = box_obstacle_zonotope('center', [3; 3; 3],...
                                     'side_lengths', [0.1; 0.1; 0.1]) ;
% obstacles{2} = box_obstacle_zonotope('center', [0.3; 0; 0.4],...
%                                      'side_lengths', [0.1; 0.8; 0.05]) ;

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
M_min_eigenvalue = 8.29938; % matlab doesn't import these from urdf so hard code into class

%% 

if plot_while_running
    figure(1); clf; view(3); grid on;
end

%%% for world

% simple rotation
start = [0;-pi/2;0;0;0;0;0];
goal = [pi/4;-pi/2;0;0;0;0;0];

random_start_goal = false;
if random_start_goal
    % random start and goal for graph planner (needs to be incorporated
    % properly eventually)
    finding_start = true;
    while finding_start
    
        % create random q
        q_rand = randomConfiguration(robot);
    
        q_aug = [q_rand; 0; 0; 0];
    
        % generate forward kinematics of random configuration
        q_fk_res = forward_kinematics(q_aug, params.nominal.T0, params.nominal.joint_axes);
    
        % get the euler angles of the random configuration
    %     q_ee_euler_angles = tform2eul(q_fk_res);
        q_ee_euler_angles = rotm2eul(q_fk_res(1:3,1:3), 'XYZ');
    
    %     % if valid, add to list of q's
        if abs(q_ee_euler_angles(1))<0.15 & abs(q_ee_euler_angles(2))<0.15 % & abs(q_ee_euler_angles(3))<0.03
            start = q_rand;
            finding_start = false;
        end
    end
    finding_goal = true;
    while finding_goal
    
        % create random q
        q_rand = randomConfiguration(robot);
    
        q_aug = [q_rand; 0; 0; 0];
    
        % generate forward kinematics of random configuration
        q_fk_res = forward_kinematics(q_aug, params.nominal.T0, params.nominal.joint_axes);
    
        % get the euler angles of the random configuration
    %     q_ee_euler_angles = tform2eul(q_fk_res);
        q_ee_euler_angles = rotm2eul(q_fk_res(1:3,1:3), 'XYZ');
    
    %     % if valid, add to list of q's
        if abs(q_ee_euler_angles(1))<0.15 & abs(q_ee_euler_angles(2))<0.15 % & abs(q_ee_euler_angles(3))<0.03
            goal = q_rand;
            finding_goal = false;
        end
    end
end

%% Graph Planner Waypoints

% Load Configuration Graph

% Note that these need to match and in future store the configuration and
% collision checking points in each node. 

config_struct = load('PlannerGraphResult2.mat');
q_valid_list = config_struct.q_valid_list;

graph_struct = load('QGraph700kPartial_AllInfo.mat');
Q_Graph = graph_struct.Q_Graph;

% Reduce Graph to Largest Connected Component

[bin, binsize] = conncomp(Q_Graph,'Type','weak');

idx = (binsize(bin) == max(binsize));

Q_Graph_Reduced = subgraph(Q_Graph,idx);
q_valid_list_reduced = q_valid_list(:,idx);

% Generate Path

% collision checking to remove nodes in collision

% find closest nodes to start and goal? assign as first and last waypoints?
% or just add start and goal as nodes, assuming they are collision free?
% but they wouldn't be connected.

% iterate through the whole graph and check for closest node? or check for
% node within a certain range?
start_dist = inf;
start_node = NaN;
goal_dist = inf;
goal_node = NaN;

for i = 1:numnodes(Q_Graph_Reduced) % length(q_valid_list_reduced)

    if norm(q_valid_list_reduced(:,i)-start) < start_dist
        start_dist = norm(q_valid_list_reduced(:,i)-start);
        start_node = i;
    end
    if norm(q_valid_list_reduced(:,i)-goal) < goal_dist
        goal_dist = norm(q_valid_list_reduced(:,i)-goal);
        goal_node = i;
    end

end
if goal_node == start_node
    fprintf('Start and Goal Nodes are the Same!!!')
end
if start_dist > 7*k_range 
    fprintf('Start Node is Far Away!!!')
    start_dist
end
if goal_dist > 7*k_range
    fprintf('Goal Node is Far Away!!!')
    goal_dist
end

% if nodes are far away, could straight line plan to them?

% generate path
[Path, dist, waypoint_nodes] = shortestpath(Q_Graph_Reduced,start_node,goal_node);
% NOTE THIS GENERATES NODES THAT CANNOT BE COMPARED TO ORIGINAL Q_GRAPH
% WITHOUT WORK
% need to add q configurations to Q_Graph and then do subgraph

% % get actual waypoint configurations
% counter = 1;
% iter = true;
% while iter == true
% 
%     % use idx to find correct waypoints
% 
% 
% end
% % waypoints = q_valid_list(:,waypoint_nodes);

Q_Nodes = table2array(Q_Graph_Reduced.Nodes)';

waypoints = Q_Nodes(:,Path);

%%

% run loop
W = kinova_grasp_world_static('create_random_obstacles_flag', false, 'goal_radius', goal_radius, 'N_obstacles',length(obstacles),'dimension',dimension,'workspace_goal_check', 0, 'verbose',verbosity, 'start', start, 'goal', goal, 'obstacles', obstacles, 'goal_type', goal_type, 'grasp_constraint_flag', true,'ik_start_goal_flag', true,'u_s', u_s, 'surf_rad', surf_rad) ;

tic;

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
                 'transmision_inertia', transmision_inertia,...
                 't_total', DURATION);

% LLC
if use_robust_input
    A.LLC = uarmtd_robust_CBF_LLC('verbose', verbosity, ...
                                  'use_true_params_for_robust', false, ...
                                  'V_max', LLC_V_max, ...
                                  'alpha_constant', alpha_constant, ...
                                  'Kr', Kr, ...
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
                   'use_cuda', use_cuda_flag,...
                   'plot_HLP_flag', true, ...
                   'lookahead_distance', 0.1, ...
                   'DURATION', DURATION) ; % 't_move_temp', t_move't_plan', t_plan,...'t_stop', t_stop % _wrapper

if if_use_RRT
    P.HLP = arm_end_effector_flat_RRT_star_HLP('plot_waypoint_flag',plot_waypoint_flag,...
                                               'plot_waypoint_arm_flag',plot_waypoint_arm_flag,...
                                               'grow_tree_mode',HLP_grow_tree_mode,...
                                               'buffer',0.1) ;
end

if if_use_graph_planner
    P.HLP = robot_arm_graph_planner_HLP('waypoints',waypoints);
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
                    'stop_sim_when_ultimate_bound_exceeded', false) ; 

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

% Calculating the Acceleration

joint_angles = A.state(A.joint_state_indices,:);
joint_angular_velocity = A.state(A.joint_speed_indices,:);

qdd_post = zeros(7,length(A.time));
% calculating the acceleration in post to compare with what is stored
for i = 2:length(A.time)
    [M, C, g] = A.calculate_dynamics(joint_angles(:,i), joint_angular_velocity(:,i), A.params.true);

    for j = 1:A.n_inputs
        M(j,j) = M(j,j) + A.transmision_inertia(j);
    end
    % can I call u=A.LLC.get_control_inputs() here with the P.info?
    
    qdd_post(:,i) = M\(A.input(:,i)-C*joint_angular_velocity(:,i)-g);
%         qdd_post(:,i) = M\(A.input(:,i)-C*joint_angular_velocity(:,i)-g);

end

figure(1001)
hold on
title('Acceleration')
% plot(A.time,qdd_post,'o')
plot(A.time,A.reference_acceleration)

% Calling RNEA

for i = 1:length(A.time)

    % clear relevant variables
    clear tau f n

    % call rnea
    [tau, f, n] = rnea(joint_angles(:,i)',joint_angular_velocity(:,i)',joint_angular_velocity(:,i)',qdd_post(:,i)',true,A.params.true); % A.acceleration(:,i)'
%     [tau, f, n] = rnea(joint_angles(:,i)',joint_angular_velocity(:,i)',joint_angular_velocity(:,i)',accel_short(:,i)',true,A.params.true); % A.acceleration(:,i)'

    % store rnea results
    Tau{i} = tau;
    F{i} = f;
    N{i} = n;

    % store the contact forces
    force(:,i) = f(:,10);
    % store the contact moments
    moment(:,i) = n(:,10);

end

% Plotting the forces

% if plot_force

    figure(3)
    % plot the x-axis force (in the tray frame)
    subplot(1,3,1)
    hold on
    plot(A.time(1:end), force(1,:), 'k')
    xlabel('time (s)')
    ylabel('x-axis Force (N)')
    axis('square')
    grid on
    % plot the y-axis force (in the tray frame)
    subplot(1,3,2)
    hold on
    plot(A.time(1:end), force(2,:), 'k')
    xlabel('time (s)')
    ylabel('y-axis Force (N)')
    axis('square')
    grid on
    % plot the z-axis force (in the tray frame)
    subplot(1,3,3)
    hold on
    plot(A.time(1:end), force(3,:), 'k')
    xlabel('time (s)')
    ylabel('z-axis Force (N)')
    axis('square')
    grid on
    
% end

% Calculating Constraints

% separation constraint
sep = -1*force(3,:);

% slipping constraint
slip = sqrt(force(1,:).^2+force(2,:).^2) - u_s.*abs(force(3,:));
slip2 = force(1,:).^2+force(2,:).^2 - u_s^2.*force(3,:).^2;

for i = 1:length(A.time)
    % tipping constraint the normal way
    ZMP_top = cross([0;0;1],moment(:,i)); % normal vector should come first
    ZMP_bottom = dot([0;0;1],force(:,i));
    ZMP(:,i) = ZMP_top/ZMP_bottom;
    ZMP_rad(i) = sqrt(ZMP(1,i)^2+ZMP(2,i)^2);
    tip(i) = ZMP_rad(i) - surf_rad;
    % tipping constraint the PZ way
    ZMP_top2 = cross([0;0;1],moment(:,i));
    ZMP_bottom2 = dot([0;0;1],force(:,i));
    tip2(i) = ZMP_top(1)^2 + ZMP_top(2)^2 - ZMP_bottom^2*(surf_rad)^2;
end

% plotting constraints

% if plot_constraint
    
    figure(4)
    % plot the separation constraint
    subplot(1,3,1)
    hold on
    plot(A.time(1:end),sep, 'k')
    xlabel('time (s)')
    ylabel('Separation Constraint')
    axis('square')
    grid on
    % plot the slipping constraint
    subplot(1,3,2)
    hold on
    plot(A.time(1:end),slip, 'k')
    plot(A.time(1:end),slip2, 'c')
    xlabel('time (s)')
    ylabel('Slipping Constraint')
    axis('square')
    grid on
    % plot the tipping constraint
    subplot(1,3,3)
    hold on
    plot(A.time(1:end),tip, 'k')
    plot(A.time(1:end),tip2, 'c')
    xlabel('time (s)')
    ylabel('Tipping Constraint')
    axis('square')
    grid on

% end

%% Mean Planning Time
plan_time = [];
for i=1:length(P.info.planning_time)
    plan_time = [plan_time P.info.planning_time{i}];
end
mean_plan_time = mean(plan_time)
max_plan_time = max(plan_time)
iterations = summary.total_iterations
num_brakes = sum(summary.stop_check)

%% numerically calculating acceleration

% numerical_accel = gradient(A.time)./gradient(A.state(A.joint_speed_indices,:));
% 
% figure()
% plot(A.time,numerical_accel)