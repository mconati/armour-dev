%% description
% This script iterates through a list of presaved random worlds and runs
% the aRmTD planner on them. It then saves information on how well each the
% planner performed in each trial.
%
% Authors: Bohao Zhang (adapted from Patrick Holmes code)
% Created 25 November 2019
% Edited 16 January 2020
% Edited: 02 June 2022 to work with updated UARMTD code
% Edited: 25 September 2022 to work with updated UARMTD code for kinova
% Edited: 08 November 2022 (by Patrick) to update the "hard scenarios"

close all; clear; clc; figure(1); clf; view(3); grid on;

%% user parameters

use_robust_input = true;

% goal_type = 'end_effector_location';
goal_type = 'configuration';
goal_radius = pi/30;
dimension = 3 ;
verbosity = 10;
% nLinks = 3 ;

%%% for planner
traj_type = 'orig';
allow_replan_errors = true ;
first_iter_pause_flag = false;
use_q_plan_for_cost = false; % otherwise use q_stop (q at final time)
input_constraints_flag = true;

%%% for agent
agent_urdf = 'gen3.urdf';

add_uncertainty_to = 'all'; % choose 'all', 'link', or 'none'
links_with_uncertainty = {}; % if add_uncertainty_to = 'link', specify links here.
% links_with_uncertainty = {'dumbbell_link'}; % if add_uncertainty_to = 'link', specify links here.
uncertain_mass_range = [0.97, 1.03];

agent_move_mode = 'integrator' ; % pick 'direct' or 'integrator'
use_CAD_flag = false;
add_measurement_noise_ = false;
measurement_noise_size_ = 0;

%%% for LLC
LLC_V_max = 5e-5;
use_true_params_for_robust = false;
if_use_mex_controller = true;

no_obstacles = false;

if_use_RRT = false;
% create_random_obstacles_flag = false ;
% verbosity = 6 ;
% actual_t_plan = 10 ;
% simulated_t_plan = 0.5 ;
% HLP_timeout = 2 ; 
HLP_grow_tree_mode = 'new' ;
% plot_while_sampling_flag = false ;
% make_new_graph_every_iteration = false ;
% plot_HLP_flag = true ; % for planner
plot_waypoint_flag = true ; % for HLP
plot_waypoint_arm_flag  = true ; % for HLP
% lookahead_distance = 0.2 ;
% use_end_effector_for_cost_flag = true ;
% plot_CAD_flag = false ; % plot the faaaaancy arm :)

% plotting
plot_while_running = false ;
% agent_camera_distance = 3 ; % default is 3
% agent_camera_position = [-3;0;1] ; % default is [-3;0;1.5]
% plot_agent_view = 'behind' ; % none, behind, above, or onboard
% plot_zonotopes = true ;

% simulation
% start_idx = 5 ;
% end_idx = 500 ;
% verbosity = 10 ;
% max_sim_time = 6000 ;
max_sim_time = 172800 ; % 48 hours
max_sim_iter = 600 ;
stop_threshold = 4 ; % number of failed iterations before exiting
% first_iter_pause_flag = false;

% file handling
save_file_header = 'trial_' ;
file_location = '../simulator_files/testing/saved_worlds/20221023_armtd' ;
if use_robust_input
    file_location = [file_location, '_robust'];
else
    file_location = [file_location, '_nominal_only'];
end
if ~exist(file_location, 'dir')
    mkdir(file_location);
end

% world file
world_file_header = 'scene';
world_file_folder = '../simulator_files/testing/saved_worlds/20221023/';
world_file_location = sprintf('%s*%s*', world_file_folder, world_file_header);
world_file_list = dir(world_file_location);


%% robot params:
robot = importrobot(agent_urdf);
robot.DataFormat = 'col';
robot.Gravity = [0 0 -9.81];
model = create_model_from_urdf(agent_urdf);
% model = rmfield(model, 'transmissionInertia');
model = rmfield(model, 'damping');
model = rmfield(model, 'friction');
params = load_robot_params(robot, ...
                           'add_uncertainty_to', add_uncertainty_to, ...
                           'links_with_uncertainty', links_with_uncertainty,...
                           'uncertain_mass_range', uncertain_mass_range);
joint_speed_limits = [-1.3963, -1.3963, -1.3963, -1.3963, -1.2218, -1.2218, -1.2218;
                       1.3963,  1.3963,  1.3963,  1.3963,  1.2218,  1.2218,  1.2218]; % matlab doesn't import these from urdf
joint_input_limits = [-56.7, -56.7, -56.7, -56.7, -29.4, -29.4, -29.4;
                       56.7,  56.7,  56.7,  56.7,  29.4,  29.4,  29.4]; % matlab doesn't import these from urdf
M_min_eigenvalue = 5.095620491878957;

use_cuda_flag = true;

%% automated from here
% run loop
tic
for idx = 1:length(world_file_list)
    clc; 
    fprintf("THIS IS WORLD %d\n\n", idx);

    % read world CSV to get start and goal, populate obstacles:
    world_filename = world_file_list(idx).name;
    [start, goal, obstacles] = load_saved_world([world_file_folder world_filename]);
%     disp('**********HACK********** pad start/goal with zeros for orig ARMTD worlds!!')
%     start = [start; 0];
%     goal = [goal; 0];
    if no_obstacles
        obstacles = {};
        disp('Obstacles disabled!!!');
        include_base_obstacle = false;
    else
        include_base_obstacle = true;
    end
    
%     figure(1); clf; view(3); grid on;
    
    W = kinova_world_static('create_random_obstacles_flag', false, 'include_base_obstacle', include_base_obstacle, 'goal_radius', goal_radius, 'N_obstacles',length(obstacles),'dimension',dimension,'workspace_goal_check', 0,...
    'verbose',verbosity, 'start', start, 'goal', goal, 'obstacles', obstacles, 'goal_type', goal_type) ;
    
%     W = fetch_base_world_static('create_random_obstacles_flag', false, 'include_base_obstacle', false, 'goal_radius', goal_radius, 'N_obstacles',length(obstacles),'dimension',dimension,'workspace_goal_check', 0,...
%     'verbose',verbosity, 'start', start, 'goal', goal, 'obstacles', obstacles, 'goal_type', goal_type) ;

    % create arm agent
    A = uarmtd_agent(robot, model, params,...
                     'verbose', verbosity,...
                     'animation_set_axes_flag', 0,... 
                     'animation_set_view_flag', 0,...
                     'move_mode', agent_move_mode,...
                     'use_CAD_flag', use_CAD_flag,...
                     'joint_speed_limits', joint_speed_limits, ...
                     'joint_input_limits', joint_input_limits, ...
                     'add_measurement_noise_', add_measurement_noise_, ...
                     'measurement_noise_size_', measurement_noise_size_,...
                     'M_min_eigenvalue', M_min_eigenvalue);

    % LLC
    if use_robust_input
        A.LLC = uarmtd_robust_CBF_LLC('verbose', verbosity, ...
                                  'use_true_params_for_robust', use_true_params_for_robust, ...
                                  'V_max', LLC_V_max, ...
                                  'if_use_mex_controller', if_use_mex_controller);
    else
        A.LLC = uarmtd_nominal_passivity_LLC('verbose', verbosity);
    end
%     A.LLC = uarmtd_robust_CBF_MEX_LLC('verbose', verbosity, ...
%                                   'use_true_params_for_robust', use_true_params_for_robust);
    A.LLC.setup(A);
    
    P = uarmtd_planner('verbose', verbosity, ...
    'first_iter_pause_flag', first_iter_pause_flag, ...
    'use_q_plan_for_cost', use_q_plan_for_cost, ...
    'input_constraints_flag', input_constraints_flag, ...
    'use_robust_input', use_robust_input, ...
    'traj_type', traj_type, ...
    'use_cuda', use_cuda_flag) ;

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
    % W.start = zeros(6, 1); % put in "home" config
    A.state(A.joint_state_indices) = W.start ;
    
    % create simulator
    S = simulator_armtd(A,W,P, 'verbose', verbosity, ...
        'stop_threshold', stop_threshold, ...
        'plot_while_running', plot_while_running,...
        'allow_replan_errors',allow_replan_errors,...
        'max_sim_time',max_sim_time,...
        'max_sim_iterations',max_sim_iter,...
        'stop_sim_when_ultimate_bound_exceeded', use_robust_input) ; 
    
    % %% plotting
    if plot_while_running
        figure(1) ; clf ; axis equal ; xlim([-1 1]); ylim([-1 1]); zlim([0 2]); grid on; hold on ;

        if dimension == 3
            view(3)
        end
        
        plot(A)
        plot(W)
    end
    
%     animate(A)
    
    % run simulation
    summary = S.run() ;
    
    % save summary
%     filename = [file_location,save_file_header,num2str(idx,'%04.f'),'.mat'] ;
    filename = [file_location,'/',save_file_header,world_filename(1:end-4),'.mat'] ;
%     save(filename, 'world_filename', 'summary')
    parsave(filename, world_filename, summary)
    
%     toc
end

function [] = parsave(filename, world_filename, summary)
    save(filename, 'world_filename', 'summary');
end