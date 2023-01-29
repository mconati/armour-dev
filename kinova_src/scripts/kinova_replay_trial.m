%% user parameters
clear;
filename = 'trial_scene_021_.mat';

verbosity = 0 ;
dimension = 3 ;

plot_start_and_end_config_only = true; % otherwise, animate trial.

%% automated from here
load(filename)

u_s = W.u_s;

agent_info = summary.agent_info ;
bounds = summary.bounds ;
obstacles = summary.obstacles ;
planner_info = summary.planner_info ;
start = summary.start ;
goal = summary.goal ;

% agent just for visualizing, parameters may differ
agent_urdf = 'Kinova_Grasp_URDF.urdf';
robot = importrobot(agent_urdf);
robot.DataFormat = 'col';
robot.Gravity = [0 0 -9.81];
params = load_robot_params(robot);

% create arm agent
% A = uarmtd_agent(robot, params,...
%              'verbose', verbosity,...
%              'animation_set_axes_flag', 0,... 
%              'animation_set_view_flag', 0);

% create world
goal_type = 'configuration';
goal_radius = pi/30;
% W = kinova_grasp_world_static('create_random_obstacles_flag', false, 'include_base_obstacle', true, 'goal_radius', goal_radius, 'N_obstacles',length(obstacles),'dimension',dimension,'workspace_goal_check', 0,...
%                             'verbose',verbosity, 'start', start, 'goal', goal, 'obstacles', obstacles, 'goal_type', goal_type,...
%                             'grasp_constraint_flag', true,'ik_start_goal_flag', true,'u_s', u_s, 'surf_rad', surf_rad) ;

% fill in agent state
A.time = agent_info.time ;
A.state = agent_info.state ;

% set up world using arm
W.setup(agent_info)
clear_plot_data(W);
    
%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(W) ;

if dimension == 3
    view(3) ;
end

if plot_start_and_end_config_only
    plot_at_time(A, 0);
    disp('Press a key to plot final config.');
    pause();
    plot(A) ;
else
    animate(A) ;
end