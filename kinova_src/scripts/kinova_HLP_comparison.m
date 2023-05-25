%% HLP Comparison

clear all;
close all;
clc;

%% Load Data

filenames = {'Simulation_w_GraphOnly.mat',...
    'Simulation_w_SLPandGraph.mat',...
    'Simulation_w_GraphOnly_LA_0p18.mat',...
    'Simulation_w_SLPandGraph_LA_0p18.mat',...
    'Simulation_w_SLPandGraphandNoise_LA_0p18.mat',...
    'Simulation_w_SLPandGraph_v2.mat',...
    'Simulation_v2.mat'};
labels = {'Graph Planner Only',...
    'Graph Planner + SLP',...
    'Graph Planner + Lookahead',...
    'Graph Planner + SLP + Lookahead',...
    'Graph Planner + SLP + Noise + Lookahead',...
    'Graph Planner + SLP + Small LA + Large WD',...
    'GP + SLP + Small LA + Large WD + No Obs'};

for i = 1:length(filenames)
    data{i} = load(filenames{i});
end
% data{1} = load('Simulation_w_GraphOnly.mat');
% data{2} = load('Simulation_w_SLPandGraph.mat');
% data{3} = load('Simulation_w_GraphOnly_LA_0p18.mat');
% data{4} = load('Simulation_w_SLPandGraph_LA_0p18.mat');
% data{5} = load('Simulation_w_SLPandGraphandNoise_LA_0p18.mat');

%% Plotting

%% Plot Trajectories

for i = 1:length(filenames)
    figure(i)
    subplot(3,1,1)
    hold on
    plot(data{i}.A.time,data{i}.A.state(data{i}.A.joint_state_indices,:))
    title(labels{i})
    subplot(3,1,2)
    hold on
    plot(data{i}.A.time,data{i}.A.state(data{i}.A.joint_speed_indices,:))
    subplot(3,1,3)
    hold on
    plot(data{i}.A.time,data{i}.A.reference_acceleration)

end

%% Plot Robot Animation

% clear all;

for i = 7:length(filenames)

    load(filenames{i})
    verbosity = 0 ;
    dimension = 3 ;
    plot_start_and_end_config_only = false; % otherwise, animate trial.
    
    A.animation_gif_filename = 'Test1.gif';
    
    u_s = W.u_s;
    agent_info = summary.agent_info ;
    bounds = summary.bounds ;
    obstacles = summary.obstacles ;
    planner_info = summary.planner_info ;
    start = summary.start ;
    goal = summary.goal ;
    
    % agent just for visualizing, parameters may differ
    agent_urdf = 'Kinova_Grasp_w_Tray.urdf';
    robot = importrobot(agent_urdf);
    robot.DataFormat = 'col';
    robot.Gravity = [0 0 -9.81];
    params = load_robot_params(robot);
    
    % create arm agent
    A = uarmtd_agent(robot, params,...
                 'verbose', verbosity,...
                 'animation_set_axes_flag', 0,... 
                 'animation_set_view_flag', 0);
    
    % create world
    goal_type = 'configuration';
    goal_radius = pi/30;
    W = kinova_grasp_world_static('create_random_obstacles_flag', false, 'include_base_obstacle', false, 'goal_radius', goal_radius, 'N_obstacles',length(obstacles),'dimension',dimension,'workspace_goal_check', 0,...
                                'verbose',verbosity, 'start', start, 'goal', goal, 'obstacles', obstacles, 'goal_type', goal_type,...
                                'grasp_constraint_flag', true,'ik_start_goal_flag', true,'u_s', u_s, 'surf_rad', surf_rad) ;
    
    % fill in agent state
    A.time = agent_info.time ;
    A.state = agent_info.state ;
    A.use_CAD_flag = true;
    n_links_and_joints = params.nominal.num_joints;
    n_states = 2*params.nominal.num_q;
    q_index = params.nominal.q_index;
    joint_state_indices = 1:2:n_states ;
    joint_speed_indices = 2:2:n_states ;
    joint_types = params.nominal.joint_types';
    joint_axes = params.nominal.joint_axes;
    link_shapes = repmat({'cuboid'}, 1, n_links_and_joints);
    [link_poly_zonotopes, link_sizes, temp_link_CAD_data] = create_pz_bounding_boxes(robot);
    A.load_CAD_arm_patch_data(temp_link_CAD_data)
    A.link_plot_edge_color = [0 0 1] ;
    A.link_plot_edge_opacity = 0 ;
    A.link_plot_face_color = [0.8 0.8 1] ;
    A.link_plot_face_opacity = 1 ;
    
    % set up world using arm
    W.setup(agent_info)
    clear_plot_data(W);
    
    % plotting
    figure(101+i) ; clf ; axis equal ; hold on ; grid on ;
    
    plot(W) ;
    
    if dimension == 3
        view(3) ;
        axis equal
    end

    P.HLP.plot_data.previous_waypoint_arm_volume = patch([1e-3 3e-3 2e-3],[1e-3 2e-3 3e-3],[1 0 0]);
    P.HLP.plot_data.next_waypoint_arm_volume = patch([1e-3 3e-3 2e-3],[1e-3 2e-3 3e-3],[1 0 0]);
    
    if plot_start_and_end_config_only
        plot_at_time(A, 0);
        grid off
        axis off
        disp('Press a key to plot final config.');
        pause();
        plot(A) ;
    else
        view(85,10)
        zoom(2)
        xlim([-0.5 1.5])
        ylim([-1.3 0.9])
        fprintf('Paused. Please resize if desired.')
        pause()
    
        animate_w_waypoint(A,P,true);
    end
end