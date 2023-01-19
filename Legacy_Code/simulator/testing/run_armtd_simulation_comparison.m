%% description
% This script iterates through a list of presaved random worlds and runs
% the aRmTD planner on them. It then saves information on how well each the
% planner performed in each trial.
%
% Authors: Patrick Holmes (adapted from Shreyas Kousik code)
% Created 25 November 2019
% Edited 16 January 2020

clear; clc;

%% user parameters
% goal_type = 'end_effector_location';
goal_type = 'configuration';
goal_radius = pi/30;
dimension = 3 ;
nLinks = 3 ;
allow_replan_errors = true ;
t_plan = 0.5 ;
time_discretization = 0.01 ;
T = 1 ;
use_cuda_flag = true;
use_new_cuda_version = true;
agent_move_mode = 'integrator';

create_random_obstacles_flag = false ;
verbosity = 6 ;
actual_t_plan = 10 ;
simulated_t_plan = 0.5 ;
HLP_timeout = 2 ; 
HLP_grow_tree_mode = 'new' ;
plot_while_sampling_flag = false ;
make_new_graph_every_iteration = false ;
plot_HLP_flag = true ; % for planner
plot_waypoint_flag = true ; % for HLP
plot_waypoint_arm_flag  = true ; % for HLP
lookahead_distance = 0.2 ;
use_end_effector_for_cost_flag = true ;
plot_CAD_flag = false ; % plot the faaaaancy arm :)

% plotting
plot_while_running = false ;
% agent_camera_distance = 3 ; % default is 3
% agent_camera_position = [-3;0;1] ; % default is [-3;0;1.5]
% plot_agent_view = 'behind' ; % none, behind, above, or onboard
% plot_zonotopes = true ;

% simulation
% start_idx = 5 ;
% end_idx = 500 ;
verbosity = 10 ;
max_sim_time = 6000 ;
max_sim_iter = 1000 ;
first_iter_pause_flag = false;

% requires a ros network be created
rosinit_fetch() ; 

% file handling
save_file_header = 'trial_' ;
file_location = '../armtd-dev/simulator_files/testing/trial_data/20211122_new_robust_input_no_obstacles_20220127' ;
if ~exist(file_location, 'dir')
    mkdir(file_location);
end

% world file
world_file_header = 'scene';
% world_file_folder = '../armtd-dev/simulator_files/testing/saved_worlds/20200116_scenarios/';
world_file_folder = '../armtd-dev/simulator_files/testing/saved_worlds/20211122/';
world_file_location = sprintf('%s*%s*', world_file_folder, world_file_header);
world_file_list = dir(world_file_location);


%% automated from here
% run loop
tic
previously_failed_trials = [89] ; 
for idx = previously_failed_trials %1:length(world_file_list)
    
    % read world CSV to get start and goal, populate obstacles:
    world_filename = world_file_list(idx).name;
    [start, goal, obstacles] = load_saved_world([world_file_folder world_filename]);
    if no_obstacles
        obstacles = {};
        disp('Obstacles disabled!!!');
        include_base_obstacle = false;
    else
        include_base_obstacle = true;
    end
    
%     figure(1); clf; view(3); grid on;

    % Create CAIK object
    fetch_CAIK = CAIK_fetch_6dof() ;

    fetch_CAIK.options_fmincon = optimoptions('fmincon', 'Display', 'off', ...
        'SpecifyObjectiveGradient', true, 'CheckGradients', false, ...
        'SpecifyConstraintGradient', true) ;

    fetch_CAIK.add_obstacle_half_plane_many(obstacles, 0.085) ; % 0.07m buffer length for obstacles
    fetch_CAIK.q_now = start ;

    J_des = [1; 0; 0.8];
    R_des = [0, 0, 1;...
             0, 1, 0;...
            -1, 0, 0];

    % Create goal configuration using CAIK
    [q_opt_0, ~] = fetch_CAIK.inverse_kinematics_obs_aware_q(J_des, R_des) ; % this is just done so the obstacles halfplanes are set up.
    fetch_CAIK.q_now = q_opt_0 ;
    
    % Create intermediate high-level waypoints using biRRT
    PRM_options = {} ;
    PRM_options.t_timeout_PRM = 10 ;        % [s]
    PRM_options.n_poses = 250000 ;          % things can go, like, pretty slow :(
    PRM_options.d_connect = 1 ;             % connectivity distance IN CONFIGURATION SPACE oolala
    PRM_options.d_discretize_for_collision_check = 0.3 ;
    PRM_options.n_neighbors = 10 ;          % min number of neighbors to connect each pose to

    PRM = PRM_planner(start, goal, obstacles, fetch_CAIK, PRM_options) ;
    PRM.generate_path() ; % generate path once (static obstacles, no further PRM computations)
    manual_waypoints = PRM.path ;    


    W = fetch_base_world_static('create_random_obstacles_flag', false, 'include_base_obstacle', true, 'goal_radius', goal_radius, 'N_obstacles',length(obstacles),'dimension',dimension,'workspace_goal_check', 0,...
    'verbose',verbosity, 'start', start, 'goal', goal, 'obstacles', obstacles, 'goal_type', goal_type) ;
    
%     W = fetch_base_world_static('create_random_obstacles_flag', false, 'include_base_obstacle', false, 'goal_radius', goal_radius, 'N_obstacles',length(obstacles),'dimension',dimension,'workspace_goal_check', 0,...
%     'verbose',verbosity, 'start', start, 'goal', goal, 'obstacles', obstacles, 'goal_type', goal_type) ;

    % create arm agent
%     A = robot_arm_3D_fetch_rnea('verbose', verbosity, 'animation_set_axes_flag', 0, 'animation_set_view_flag', 0, 'move_mode', agent_move_mode);
    A = robot_arm_3D_fetch_rnea_new_robust('verbose', verbosity, 'animation_set_axes_flag', 0, 'animation_set_view_flag', 0, 'move_mode', agent_move_mode);

    
    % can adjust LLC gains here
    A.LLC.K_p = 1*A.LLC.K_p;
    A.LLC.K_i = 1*A.LLC.K_i;
    A.LLC.K_d = 1*A.LLC.K_d;
%     A.joint_input_limits = 1*A.joint_input_limits;
    
    % options for FRS
    FRS_options = struct();
    FRS_options.t_plan = t_plan;
    FRS_options.origin_shift = A.joint_locations(1:3, 1);
    FRS_options.T = T;
    FRS_options.buffer_dist = A.buffer_dist;
    FRS_options.combs = generate_combinations_upto(200);
    FRS_options.maxcombs = 200;
    
    % create planner (OLD)
%     P = robot_arm_rotatotope_RTD_planner_3D_fetch(FRS_options,...
%         'verbose', verbosity, 't_plan', actual_t_plan,...
%         't_move',simulated_t_plan,...
%         'time_discretization', time_discretization,...
%         'first_iter_pause_flag',first_iter_pause_flag,...
%         'plot_HLP_flag',plot_HLP_flag,...
%         'lookahead_distance',lookahead_distance,...
%         'use_end_effector_for_cost_flag',use_end_effector_for_cost_flag,...
%         'use_cuda_flag', use_cuda_flag) ;

    % create planner (NEW 20211117)
%     P = robot_rotatotope_RTD_planner('fetch_rotatotope_FRS', FRS_options,...
%         'verbose', verbosity, 't_plan', actual_t_plan,...
%         't_move',simulated_t_plan,...
%         'time_discretization', time_discretization,...
%         'first_iter_pause_flag',first_iter_pause_flag,...
%         'plot_HLP_flag',plot_HLP_flag,...
%         'lookahead_distance',lookahead_distance,...
%         'use_end_effector_for_cost_flag',use_end_effector_for_cost_flag,...
%         'use_cuda_flag', use_cuda_flag,...
%         'use_new_cuda_version', use_new_cuda_version) ;
    
%     P = robot_RTD_noFRS_planner('verbose', verbosity, 't_plan', actual_t_plan,...
%         't_move',simulated_t_plan,...
%         'time_discretization', time_discretization,...
%         'first_iter_pause_flag',first_iter_pause_flag,...
%         'plot_HLP_flag',plot_HLP_flag,...
%         'lookahead_distance',lookahead_distance,...
%         'use_end_effector_for_cost_flag',use_end_effector_for_cost_flag,...
%         'use_cuda_flag', use_cuda_flag,...
%         'use_new_cuda_version', use_new_cuda_version) ;
    
    P = robot_RTD_noFRS_planner('verbose', verbosity, 't_plan', actual_t_plan,...
        't_move',simulated_t_plan,...
        'time_discretization', time_discretization,...
        'first_iter_pause_flag',first_iter_pause_flag,...
        'plot_HLP_flag',plot_HLP_flag) ;

    
%     P.HLP = arm_end_effector_RRT_star_HLP('plot_waypoint_flag',plot_waypoint_flag,...
%         'plot_waypoint_arm_flag',plot_waypoint_arm_flag,...
%         'grow_tree_mode',HLP_grow_tree_mode,...
%         'buffer',0.1) ;
    
    P = robot_PRM_RTD_planner('fetch_rotatotope_FRS', FRS_options,...
        'verbose', verbosity,...
        't_plan', t_plan,...
        'time_discretization', time_discretization,...
        'plot_HLP_flag',plot_HLP_flag,...
        'use_cuda_flag', use_cuda_flag,...
        'use_new_cuda_version', use_new_cuda_version,...
        'first_iter_pause_flag', first_iter_pause_flag) ;

    P.HLP = arm_iterative_manual_waypoint_HLP('manual_waypoints',manual_waypoints) ; 


    % set up world using arm
    I = A.get_agent_info ;
    W.setup(I)
    
    % place arm at starting configuration
    % W.start = zeros(6, 1); % put in "home" config
    A.state(A.joint_state_indices) = W.start ;
    
    % create simulator
    S = simulator_armtd(A,W,P, 'plot_while_running', plot_while_running, 'allow_replan_errors',allow_replan_errors,'max_sim_time',max_sim_time,'max_sim_iterations',max_sim_iter) ; 
    
    % %% plotting
    if plot_while_running
        figure(1) ; clf ; axis equal ; hold on ; grid on
        
        plot(A)
        plot(W)
        
        if dimension == 3
            view(3)
        end
    end
    
%     animate(A)
    
    % run simulation
    summary = S.run() ;
    
    % save summary
%     filename = [file_location,save_file_header,num2str(idx,'%04.f'),'.mat'] ;
    filename = [file_location,'/',save_file_header,world_filename(1:end-4),'.mat'] ;
    save(filename, 'world_filename', 'summary')
    
    toc
end

