%% user parameters
<<<<<<< HEAD
clear;
% close(2,3)

% filename = 'trial_scene_059_.mat';
filename = ['trial_scene_010_098.csv.mat'];
=======
close all; clear; clc;

file_location = '../results/random/' ;
scene_idx = 44;
>>>>>>> pybind

verbosity = 0 ;
dimension = 3 ;

plot_start_and_end_config_only = true; % otherwise, animate trial.

%% automated from here
<<<<<<< HEAD
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

%% Output some key stats

goal_check = summary.goal_check
num_iter = summary.total_iterations
    
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
    animate(A);
end

figure()
plot(A.time,A.reference_acceleration)
title('Joint Acceleration')
figure()
plot(A.time,A.state(A.joint_speed_indices,:))
title('Joint Speeds')
=======
load(['agent_info_',num2str(scene_idx),'.mat']);
    
%% plotting
figure(1) ; clf ; axis equal ; xlim([-1 1]); ylim([-1 1]); zlim([0 1.5]); grid on; hold on ; view(3) ;

A.animation_gif_filename = [file_location,'/','animation_',num2str(scene_idx),'.gif'];
frame_rate = A.animation_time_discretization / A.animation_playback_rate ;
start_gif = true;
W.plot();
A.plot_at_time(A.time(1));
breakpoint = 1;
for i = 1:10:(length(A.time)-1)
    A.plot_at_time(A.time(i));
   
    r_id = floor(i / 50) + 1;
    for j = 1:7
        hd{j} = trisurf(convhulln(P.info.sliced_FO_zono{r_id}{j}), ...
                        P.info.sliced_FO_zono{r_id}{j}(:,1), ...
                        P.info.sliced_FO_zono{r_id}{j}(:,2), ...
                        P.info.sliced_FO_zono{r_id}{j}(:,3), ...
                        'FaceColor',[0,0.5,0.5], ...
                        'FaceAlpha',0.1, ...
                        'EdgeColor',[0,0.5,0.5], ...
                        'EdgeAlpha',0.2);
    end

    fh = get(groot,'CurrentFigure') ;
    frame = getframe(fh) ;
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    if start_gif
        imwrite(imind,cm,A.animation_gif_filename,'gif', 'Loopcount',inf,...
                'DelayTime',frame_rate) ;
        start_gif = false ;
    else 
        imwrite(imind,cm,A.animation_gif_filename,'gif','WriteMode','append',...
                'DelayTime',frame_rate) ;
    end
    
    for j = 1:7
        delete(hd{j});
    end
end
>>>>>>> pybind
