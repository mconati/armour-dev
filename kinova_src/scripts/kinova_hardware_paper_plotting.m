%% Hardware Experiment Contact Diagram Plotting
% Zachary Brei

clear all;
close all;
clc;
fig_num = 1;

%% Description

% use 
    % debug_traj_pos
    % debug_traj_vel
    % debug_traj_accel
    % debug_traj_k
% to calculate unsliced reachsets again and also the nominal values

%% Load Hardware ROS Data

load_file = 'HardwareVideoROSData_04222023.mat';
data = load(load_file);

% contact parameters: match with hardware experiment settings
u_s = 0.6;
surf_rad = 0.058/2;

%% Extract Data

% process the raw time measurements to offset to zero
time = data.rosbag.raw_time - data.rosbag.raw_time(1);

% for nominal values
pos = data.rosbag.debug_traj_pos;
vel = data.rosbag.debug_traj_vel;
accel = data.rosbag.debug_traj_accel;
opt_k = data.rosbag.debug_traj_k;

% for reach set values
duration = data.rosbag.debug_duration; % for determining if braking maneuver occured
rs_pos = data.rosbag.traj_pos;
rs_vel = data.rosbag.traj_vel;
rs_accel = data.rosbag.traj_accel;
rs_opt_k = data.rosbag.traj_k;

%% Load Robot Parameters

% match with Hardware experiment settings
add_uncertainty_to = 'all'; % choose 'all', 'link', or 'none'
links_with_uncertainty = {}; % if add_uncertainty_to = 'link', specify links here.
uncertain_mass_range = [0.97, 1.03];

agent_urdf = 'Kinova_Grasp_w_Tray.urdf';
robot = importrobot(agent_urdf);
robot.DataFormat = 'col';
robot.Gravity = [0 0 -9.81];
params = load_robot_params(robot, ...
                           'add_uncertainty_to', add_uncertainty_to, ...
                           'links_with_uncertainty', links_with_uncertainty,...
                           'uncertain_mass_range', uncertain_mass_range);

%% Calculate Nominal RNEA Values (Time Consuming)

% actuator inertia: match with hardware experiment settings
transmision_inertia = [8.02999999999999936 11.99620246153036440 9.00254278617515169 11.58064393167063599 8.46650409179141228 8.85370693737424297 8.85873036646853151]; 

joint_angles = pos';
joint_angular_velocity = vel';

qdd_post = accel'; % zeros(7,length(A.time));
%     % calculating the acceleration in post to compare with what is stored
% for i = 1:length(A.time(1:end-1))
%     [M, C, g] = A.calculate_dynamics(joint_angles(:,i), joint_angular_velocity(:,i), params.true);
%     for j = 1:size(M,1)
%         M(j,j) = M(j,j) + transmision_inertia(j);
%     end
%     qdd_post(:,i) = M\(A.input(:,i+1)-C*joint_angular_velocity(:,i)-g);
% end

% calling rnea
for i = 1:length(pos)

    % clear relevant variables
    clear tau f n

    % call rnea
    [tau, f, n] = rnea(joint_angles(:,i)',joint_angular_velocity(:,i)',joint_angular_velocity(:,i)',qdd_post(:,i)',true,params.true); % A.acceleration(:,i)'

    % store rnea results
    Tau{i} = tau;
    F{i} = f;
    N{i} = n;

    % store the contact forces
    force(:,i) = f(:,10);
    % store the contact moments
    moment(:,i) = n(:,10);

end

%% Calculate Nominal Constraint Values

% separation constraint
sep = -1*force(3,:);

% slipping constraint
slip = sqrt(force(1,:).^2+force(2,:).^2) - u_s.*abs(force(3,:));
slip2 = force(1,:).^2+force(2,:).^2 - u_s^2.*force(3,:).^2;

for i = 1:length(pos)
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

%% Iterate Through Planning Iterations?

% need to plot all of the nominal values
% update the reach sets each time a planning iteration passes
    % extend the reach sets when a braking maneuver happens?
    % check what the duration is and plot the correct amount of reach sets
    % based on that. 

% outer loop iterates through length of pos
    % check if the time counter is greater than duration
        % if so, 
            % check next (two down?) duration value to see if next planning iteration
            % has a braking maneuver.
            % update the reach sets.
            % reset the time counter.
        % else increment time counter and continue
    % plot the nominal value

%% Calculate Unsliced Reachable Sets

% get agent info
agent_info = A.get_agent_info() ; % how to get agent info? what specifically is needed?

% get world info
% given the current state of the agent, query the world
% to get the surrounding obstacles
world_info = W.get_world_info(agent_info,P) ; % how to get world info? what specifically is needed?
% might need to hard code obstacle position and size

% initial condition
q_0 = A.state(A.joint_state_indices,traj_index);
q_dot_0 = A.state(A.joint_speed_indices,traj_index);
q_ddot_0 = A.reference_acceleration(:,traj_index);

q_des = P.HLP.waypoints(:,plan_index);

% organize input to cuda program
P.vdisp('Calling CUDA & C++ Program!',3);
% P.kinova_test_folder_path = 'C:\Users\roahmlab\Documents\GitHub\armour-dev\kinova_src';
P.kinova_test_folder_path = '/home/baiyuew/ROAHM/armour-dev/kinova_src';
% cuda_input_file = fopen([P.kinova_test_folder_path, '\kinova_simulator_interfaces\kinova_planner_realtime\buffer\armour.in'], 'w');  
cuda_input_file = fopen([P.kinova_test_folder_path, '/kinova_simulator_interfaces/kinova_planner_realtime/buffer/armour.in'], 'w');  

for ind = 1:length(q_0)
    fprintf(cuda_input_file, '%.10f ', q_0(ind));
end
fprintf(cuda_input_file, '\n');
for ind = 1:length(q_dot_0)
    fprintf(cuda_input_file, '%.10f ', q_dot_0(ind));
end
fprintf(cuda_input_file, '\n');
for ind = 1:length(q_ddot_0)
    fprintf(cuda_input_file, '%.10f ', q_ddot_0(ind));
end
fprintf(cuda_input_file, '\n');
for ind = 1:length(q_des)
    fprintf(cuda_input_file, '%.10f ', q_des(ind));
end
fprintf(cuda_input_file, '\n');
fprintf(cuda_input_file, '%d\n', max(length(world_info.obstacles), 0));
for obs_ind = 1:length(world_info.obstacles)
    temp = reshape(world_info.obstacles{obs_ind}.Z, [1,size(world_info.obstacles{obs_ind}.Z,1) * size(world_info.obstacles{obs_ind}.Z,2)]);
    for ind = 1:length(temp)
        fprintf(cuda_input_file, '%.10f ', temp(ind));
    end
    fprintf(cuda_input_file, '\n');
end

fclose(cuda_input_file);

% call cuda program in terminal
% you have to be in the proper path!
terminal_output = system('env -i bash -i -c "../kinova_simulator_interfaces/kinova_planner_realtime/rtd_force_main_v2" '); % rtd-force path
%                     terminal_output = system('env -i bash -i -c "./../kinova_simulator_interfaces/kinova_planner_realtime_original/armour_main"'); % armour path

if terminal_output == 0
    data = readmatrix([P.kinova_test_folder_path, '/kinova_simulator_interfaces/kinova_planner_realtime/buffer/armour.out'], 'FileType', 'text');
%                         data = readmatrix([P.kinova_test_folder_path, '/kinova_simulator_interfaces/kinova_planner_realtime_original/buffer/armour.out'], 'FileType', 'text');
    k_opt = data(1:end-1);
    planning_time = data(end) / 1000.0; % original data is milliseconds

    if length(k_opt) == 1
        P.vdisp('Unable to find new trajectory!',3)
        k_opt = nan;
    elseif planning_time > P.t_plan
        P.vdisp('Solver Took Too Long!',3)
        k_opt = nan;
    else
        P.vdisp('New trajectory found!',3);
        for i = 1:length(k_opt)
            fprintf('%7.6f ', k_opt(i));
        end
        fprintf('\n');
    end
else
    error('CUDA program error! Check the executable path in armour-dev/kinova_src/kinova_simulator_interfaces/uarmtd_planner');
end

if terminal_output == 0
    % read FRS information if needed
    link_frs_center = readmatrix('armour_joint_position_center.out', 'FileType', 'text');
    link_frs_generators = readmatrix('armour_joint_position_radius.out', 'FileType', 'text');
    control_input_radius = readmatrix('armour_control_input_radius.out', 'FileType', 'text');
    constraints_value = readmatrix('armour_constraints.out', 'FileType', 'text');
    contact_constraint_radii = readmatrix('armour_force_constraint_radius.out', 'FileType', 'text');
    wrench_radii = readmatrix('armour_wrench_values.out', 'FileType', 'text');
else
    k_opt = nan;
end

% parse the unsliced PZ output
[force_PZ_unsliced,moment_PZ_unsliced] = parse_unsliced_PZs_func();