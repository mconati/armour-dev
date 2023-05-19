%% Hardware Experiment Contact Diagram Plotting
% Zachary Brei

%%%%% NEEDS TO BE FIXED WITH ACCELERATION PROCESSING %%%%%%

clear all;
close all;
clc;
fig_num = 1;

%% Parameters

fontsize = 16;
linewidth = 3;
% number of time steps to skip
skip_num = 1;
% planning iteration to plot
plan_index = 1;
% start index of the planning iteration?
traj_index = 1;

% match to hardware settings
time_steps = 50; 
time_index = time_steps/2; % for plotting half a trajectory
duration = 1.75;
dt = duration/time_steps;
t_traj = linspace(0,duration,time_steps);
surf_rad = 0.058/2;
u_s = 0.6;

%% Description

% use 
    % debug_traj_pos
    % debug_traj_vel
    % debug_traj_accel
    % debug_traj_k
% to calculate unsliced reachsets again and also the nominal values

%% Color Coding

unsliced_color = 1/256*[90,200,243];
% slice_color = 1/256*[72,181,163]; % light green
slice_color = 1/256*[75,154,76]; % dark green
face_alpha = 0.3;
face_alpha_light = 0.03;
edge_alpha = 0.5;

%% Load Hardware ROS Data

load_file = 'HardwareVideo_MultipleTrials_05_17_2023_ROSData.mat';
agent_urdf = 'Kinova_Grasp_w_Tray.urdf';

% load_file = 'HardwareFailureROSData.mat';
% agent_urdf = 'Kinova_Grasp_URDF.urdf';

data = load(load_file);

% contact parameters: match with hardware experiment settings
u_s = 0.6;
surf_rad = 0.058/2;

%% Extract Data

% process the raw time measurements to offset to zero
exp_time = data.rosbag.raw_time - data.rosbag.raw_time(1);

% for nominal values
pos = data.rosbag.debug_traj_pos;
vel = data.rosbag.debug_traj_vel;
% accel = data.rosbag.debug_traj_accel; % acceleration needs to be calculated
opt_k = data.rosbag.debug_traj_k;

control_torque = data.rosbag.control_torque;

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
input = control_torque';

qdd_post = zeros(7,length(exp_time)); % accel'; % 
    % calculating the acceleration in post to compare with what is stored
exp_time_2 = exp_time(1:end-1);
parfor i = 1:length(exp_time_2)
    [M, C, g] = calculate_dynamics(joint_angles(:,i)', joint_angular_velocity(:,i)', params.true);
    for j = 1:size(M,1)
        M(j,j) = M(j,j) + transmision_inertia(j);
    end
    qdd_post(:,i) = M\(input(:,i+1)-C*joint_angular_velocity(:,i)-g);
end

%% Plotting the Acceleration

fig_num = fig_num + 1;
figure(fig_num)
hold on

plot(exp_time,qdd_post)

%%

% calling rnea
Tau = cell(1,length(pos));
F = cell(1,length(pos));
N = cell(1,length(pos));
force = zeros(3,length(pos));
moment = zeros(3,length(pos));
parfor i = 1:length(pos)

%     % clear relevant variables
%     clear tau f n

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

%% Calculate Unsliced Reachable Sets

% get agent info
% agent_info = A.get_agent_info() ; % how to get agent info? what specifically is needed?

% get world info
% given the current state of the agent, query the world
% to get the surrounding obstacles
% world_info = W.get_world_info(agent_info,P) ; % how to get world info? what specifically is needed?
% might need to hard code obstacle position and size

% ! don't need the obstacle for the unsliced reach sets ! and then can
% slice using the provided k

% initial condition
q_0 = rs_pos(1,:);
q_dot_0 = rs_vel(1,:);
q_ddot_0 = rs_accel(1,:);

% pass in the k_opt as 
q_des = rs_opt_k(1,:);

% organize input to cuda program
% P.vdisp('Calling CUDA & C++ Program!',3);
% P.kinova_test_folder_path = 'C:\Users\roahmlab\Documents\GitHub\armour-dev\kinova_src';
kinova_test_folder_path = '/home/baiyuew/ROAHM/armour-dev/kinova_src';
% cuda_input_file = fopen([P.kinova_test_folder_path, '\kinova_simulator_interfaces\kinova_planner_realtime\buffer\armour.in'], 'w');  
cuda_input_file = fopen([kinova_test_folder_path, '/kinova_simulator_interfaces/kinova_planner_realtime/buffer/armour.in'], 'w');  

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
% fprintf(cuda_input_file, '\n');
% fprintf(cuda_input_file, '%d\n', max(length(world_info.obstacles), 0));
% for obs_ind = 1:length(world_info.obstacles)
%     temp = reshape(world_info.obstacles{obs_ind}.Z, [1,size(world_info.obstacles{obs_ind}.Z,1) * size(world_info.obstacles{obs_ind}.Z,2)]);
%     for ind = 1:length(temp)
%         fprintf(cuda_input_file, '%.10f ', temp(ind));
%     end
%     fprintf(cuda_input_file, '\n');
% end

fclose(cuda_input_file);

% call cuda program in terminal
% you have to be in the proper path!
terminal_output = system('env -i bash -i -c "../kinova_simulator_interfaces/kinova_planner_realtime/rtd_force_main_v2" '); % rtd-force path
%                     terminal_output = system('env -i bash -i -c "./../kinova_simulator_interfaces/kinova_planner_realtime_original/armour_main"'); % armour path

% if terminal_output == 0
%     data = readmatrix([P.kinova_test_folder_path, '/kinova_simulator_interfaces/kinova_planner_realtime/buffer/armour.out'], 'FileType', 'text');
% %                         data = readmatrix([P.kinova_test_folder_path, '/kinova_simulator_interfaces/kinova_planner_realtime_original/buffer/armour.out'], 'FileType', 'text');
%     k_opt = data(1:end-1);
%     planning_time = data(end) / 1000.0; % original data is milliseconds
% 
%     if length(k_opt) == 1
%         P.vdisp('Unable to find new trajectory!',3)
%         k_opt = nan;
%     elseif planning_time > P.t_plan
%         P.vdisp('Solver Took Too Long!',3)
%         k_opt = nan;
%     else
%         P.vdisp('New trajectory found!',3);
%         for i = 1:length(k_opt)
%             fprintf('%7.6f ', k_opt(i));
%         end
%         fprintf('\n');
%     end
% else
%     error('CUDA program error! Check the executable path in armour-dev/kinova_src/kinova_simulator_interfaces/uarmtd_planner');
% end

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

%% Parse the Unsliced Output

% parse the unsliced PZ output
[force_PZ_unsliced,moment_PZ_unsliced] = parse_unsliced_PZs_func();

%% Process Sliced Wrench PZ

force_vertices = [];
for i = 1 % :length(sim_result.P.info.wrench_radii)

%     wrench = wrench_radii{i};
    
    % get center and radii values of sliced PZs
    force_center = wrench_radii(:,1:3);
    moment_center = wrench_radii(:,4:6);
    force_radii = wrench_radii(:,7:9);
    moment_radii = wrench_radii(:,10:12);

end

% getting upper and lower bounds for patch plotting
force_lower = force_center - force_radii;
force_upper = force_center + force_radii;

%% Calculate ZMP Overapproximation 

% unsliced
for i = 1:time_index
    if size(moment_PZ_unsliced{i}.Z,2) == 2
        % add two columns of zeros
        moment_zono = polyZonotope_ROAHM(moment_PZ_unsliced{i}.Z(:,1),[moment_PZ_unsliced{i}.Z(:,2), zeros(3,2)],[]);
    elseif size(moment_PZ_unsliced{i}.Z,2) == 3
        % add one column of zeros
        moment_zono = polyZonotope_ROAHM(moment_PZ_unsliced{i}.Z(:,1),[moment_PZ_unsliced{i}.Z(:,2), moment_PZ_unsliced{i}.Z(:,3), zeros(3,1)],[]);
    else
        moment_zono = polyZonotope_ROAHM(moment_PZ_unsliced{i}.Z(:,1),moment_PZ_unsliced{i}.Z(:,2:end),[]);
    end
    % calculating the PZ form of the constraints
    ZMP_PZ_top = cross([0;0;1],moment_zono);
    ZMP_PZ_bottom = interval(force_PZ_unsliced{i}); % *[0,0,1]; % modified
%     ZMP_PZ_bottom = interval(ZMP_PZ_bottom);
    ZMP_PZ_bottom_inf = ZMP_PZ_bottom.inf(3);
    ZMP_PZ_bottom_sup = ZMP_PZ_bottom.sup(3);
    % Note that this is not the entire overapproximated PZ because the
    % bottom has the .sup as well
    ZMP_PZ_inf = interval(ZMP_PZ_top) / ZMP_PZ_bottom_inf;
    ZMP_PZ_sup = interval(ZMP_PZ_top) / ZMP_PZ_bottom_sup;
    ZMP_unsliced_PZ{i} = convHull(ZMP_PZ_inf,ZMP_PZ_sup);
end

% sliced
for i = 1:time_index
    force_int{i} = interval(force_center(i,3)-force_radii(i,3),force_center(i,3)+force_radii(i,3));
    moment_zono = polyZonotope_ROAHM(moment_center(i,:)',eye(3).*moment_radii(i,:)');
    % calculating the PZ form of the constraints
    ZMP_PZ_top = cross([0;0;1],moment_zono);
    ZMP_PZ_bottom = force_int{i}; % *[0,0,1]; % modified
    ZMP_PZ_bottom_inf = ZMP_PZ_bottom.inf;
    ZMP_PZ_bottom_sup = ZMP_PZ_bottom.sup;
    % Note that this is not the entire overapproximated PZ because the
    % bottom has the .sup as well
    ZMP_PZ_inf = interval(ZMP_PZ_top) / ZMP_PZ_bottom_inf;
    ZMP_PZ_sup = interval(ZMP_PZ_top) / ZMP_PZ_bottom_sup;
    ZMP_PZ{i} = convHull(ZMP_PZ_inf,ZMP_PZ_sup);
end

%% Plotting Separation

fig_num = fig_num + 1;
figure(fig_num); clf; hold on;

max_sup = -1;

% plot unsliced force overapproximation
for i = 1:time_index

    force_int_unsliced = interval(force_PZ_unsliced{i});

    fz1 = patch([t_traj(i)+dt; t_traj(i)+dt; t_traj(i); t_traj(i)], [force_int_unsliced.sup(3); force_int_unsliced.inf(3); force_int_unsliced.inf(3); force_int_unsliced.sup(3)],'b');
    fz1.LineWidth = 0.01;
    fz1.FaceColor = unsliced_color;
    fz1.FaceAlpha = face_alpha;
    fz1.EdgeAlpha = edge_alpha;

    if force_int_unsliced.sup(3) > max_sup
        max_sup = force_int_unsliced.sup(3);
    end

end

% plot sliced force overapproximation
for i = 1:time_index
    fz2 = patch([t_traj(i)+dt; t_traj(i)+dt; t_traj(i); t_traj(i)], [force_upper(i,3); force_lower(i,3); force_lower(i,3); force_upper(i,3)],'b');
    fz2.LineWidth = 0.01;
    fz2.FaceColor = slice_color;
    fz2.FaceAlpha = face_alpha;
    fz2.EdgeAlpha = edge_alpha;
end

% plot the nominal values
nom_idx = find(exp_time < t_traj(time_index+1));
plot(exp_time(nom_idx), force(3,nom_idx),'-k','LineWidth',linewidth)

% formatting for plot
ylim([0 max_sup+0.1])
xlim([0 t_traj(time_index+1)])
xlabel('Time (s)')
ylabel('Force (N)')
title('Contact Joint z-Axis Force')
set(gca,'FontSize',fontsize)



%% Plotting Friction Cone Diagram

fig_num = fig_num + 1;
figure(fig_num)
hold on

% plotting friction cone boundary
max_fz = max(force(3,1:time_index)); % UPDATE THIS SHOULDN't BE TIME INDEX
min_fz = min(force(3,1:time_index));
[max_fz_x, max_fz_y] = circle(2*max_fz*u_s,0,0,0,2*pi,0.01);
[min_fz_x, min_fz_y] = circle(2*min_fz*u_s,0,0,0,2*pi,0.01);
plot(max_fz_x, max_fz_y,'-r')
plot(min_fz_x, min_fz_y,'-r')
fill([max_fz_x flip(min_fz_x)],[max_fz_y flip(min_fz_y)],'r','FaceAlpha',0.9)
fp_bound = line([max_fz_x(end) min_fz_x(1)],[max_fz_y(end) min_fz_y(1)]);
set(fp_bound,'Color','r') % ,'EdgeAlpha',0.3)

% % plotting unsliced overapproximation
% for i = 1:skip_num:time_index
% 
%     force_int_unsliced = interval(force_PZ_unsliced{i});
% 
%     patch([force_int_unsliced.inf(1),force_int_unsliced.sup(1),force_int_unsliced.sup(1),force_int_unsliced.inf(1)],[force_int_unsliced.inf(2),force_int_unsliced.inf(2),force_int_unsliced.sup(2),force_int_unsliced.sup(2)],unsliced_color,'FaceAlpha',face_alpha_light,'EdgeAlpha',edge_alpha,'EdgeColor',unsliced_color)
% 
% end
% 
% % plotting sliced overapproximation
% for i = 1:skip_num:time_index
% 
%     patch([force_lower(i,1),force_upper(i,1),force_upper(i,1),force_lower(i,1)],[force_lower(i,2),force_lower(i,2),force_upper(i,2),force_upper(i,2)],slice_color,'FaceAlpha',face_alpha_light,'EdgeAlpha',edge_alpha,'EdgeColor',slice_color)
% 
% end

% plotting nominal values
nom_idx = find(exp_time < t_traj(time_index+1));
plot(force(1,nom_idx),force(2,nom_idx),'-k')

% plot formatting
xlabel('x-axis Force (N)')
ylabel('y-axis Force (N)')
title('Friction Cone')
set(gca,'FontSize',fontsize)
axis square
axis equal

%% Plotting ZMP Diagram

fig_num = fig_num + 1;
figure(fig_num)
hold on

% scale from meter to centimeter
factor = 100; 

% plot constraint boundary
r=surf_rad*factor;
x=0;
y=0;
th = linspace(0,2*pi,500);
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
tipplot1 = plot(xunit, yunit,'-r');


% plot unsliced
for i = 1:skip_num:time_index
    zmp1 = plot(ZMP_unsliced_PZ{i}.*factor,[1,2],'Filled',true);
    zmp1.LineWidth = 0.1;
    zmp1.FaceColor = unsliced_color;
    zmp1.EdgeColor = unsliced_color;
    zmp1.EdgeAlpha = edge_alpha;
    zmp1.FaceAlpha = face_alpha_light;
end

% plot sliced
for i = 1:skip_num:time_index
    zmp2 = plot(ZMP_PZ{i}.*factor,[1,2],'Filled',true);
    zmp2.LineWidth = 0.1;
    zmp2.FaceColor = slice_color;
    zmp2.EdgeColor = slice_color;
    zmp2.EdgeAlpha = edge_alpha;
    zmp2.FaceAlpha = face_alpha_light;
end

% plot nominal
nom_idx = find(exp_time < t_traj(time_index+1));
plot(ZMP(1,nom_idx).*factor,ZMP(2,nom_idx).*factor,'-k')

% plot formatting
xlabel('x_o position (cm)')
ylabel('y_o position (cm)')
title('Zero Moment Point')
axis('square')
axis equal
grid on
set(gca,'FontSize',fontsize)

%% helper functions
function [link_poly_zonotopes, link_sizes, mesh] = create_link_poly_zonos(robot)
% takes in a robot
% creates overapproximating PZ bounding box for each body of the robot.

    for i = 1:robot.NumBodies
        bounds = zeros(3, 2);
      
        if strcmp(robot.Bodies{i}.Visuals{1}, 'primitive:box')
            col_str = robot.Bodies{i}.Collisions{1};
        
            col_str = split(col_str(11:end-1));
            
            link_sizes(:,i) = zeros(length(col_str),1);
        
            for j=1:length(col_str)
                link_sizes(j,i) = str2num(col_str{j});
            end

            c = link_sizes(:,i)/2;
            c(2) = 0;
            Grest = diag(link_sizes(:,i)/2);
    
            link_poly_zonotopes{i,1} = polyZonotope_ROAHM(c, [], Grest);
            
    
            mesh = [];
    
        elseif contains(robot.Bodies{i}.Visuals{1}, '.stl') || contains(robot.Bodies{i}.Visuals{1}, '.STL')
            stl_file = robot.Bodies{i}.Visuals{1};
            stl_file = extractAfter(stl_file, 'Mesh:');
            mesh{i, 1} = stlread(stl_file);
            
            for j = 1:3
                bounds(j, 1) = min(mesh{i}.Points(:, j));
                bounds(j, 2) = max(mesh{i}.Points(:, j));
            end
    
            c = (bounds(:, 1) + bounds(:, 2))./2;
            g = (bounds(:, 2) - bounds(:, 1))./2;
        
            centers(i,:) = c';
            
            link_poly_zonotopes{i, 1} = polyZonotope_ROAHM(c, [], [g(1), 0, 0;...
                                                                   0, g(2), 0;...
                                                                   0, 0, g(3)]);
                                                        
            link_sizes(:, i) = [2*g(1); 2*g(2); 2*g(3)];
        else
            % 10 cm cube if STL file not available
            bounds = [-0.05, 0.05; 
                      -0.05, 0.05; 
                      -0.05, 0.05];
    
    
            c = (bounds(:, 1) + bounds(:, 2))./2;
            g = (bounds(:, 2) - bounds(:, 1))./2;
        
            centers(i,:) = c';
            
            link_poly_zonotopes{i, 1} = polyZonotope_ROAHM(c, [], [g(1), 0, 0;...
                                                                   0, g(2), 0;...
                                                                   0, 0, g(3)]);
                                                        
            link_sizes(:, i) = [2*g(1); 2*g(2); 2*g(3)];
        end
    end
end

function [xlinks, ylinks] = link_points(robot, q)

    num_bodiess = robot.NumBodies;

    xlinks = zeros(1, num_bodiess);
    ylinks = zeros(1, num_bodiess);

    % robot links
    for i = 1:robot.NumBodies
        T = getTransform(robot, q, robot.Bodies{i}.Name);
        xlinks(i) = T(1,4);
        ylinks(i) = T(2,4);
    end
end

function [xee, yee] = ee_points(robot, q, ee_size)
    
    x_width = ee_size(1);
    y_width = ee_size(2);

    T = getTransform(robot, q, robot.Bodies{end}.Name);
%     base_point = T(:,4);

    p1 = T*([0; y_width/2; 0; 1]);
    p2 = T*([x_width; y_width/2; 0; 1]);
    p3 = T*([0; -y_width/2; 0; 1]);
    p4 = T*([x_width; -y_width/2; 0; 1]);

    xee = [p1(1), p2(1), p3(1), p4(1)];
    yee = [p1(2), p2(2), p3(2), p4(2)];
end

function [q_des, qd_des, qdd_des] = desired_trajectory(P, q_0, q_dot_0, q_ddot_0, t, k)
    % at a given time t and traj. param k value, return
    % the desired position, velocity, and acceleration.
    switch P.traj_type
    case 'orig'
        t_plan = P.t_plan;
        t_stop = P.t_stop;
        k_scaled = P.jrs_info.c_k_orig + P.jrs_info.g_k_orig.*k;
        
        if ~isnan(k)
            if t <= t_plan
                % compute first half of trajectory
                q_des = q_0 + q_dot_0*t + (1/2)*k_scaled*t^2;
                qd_des = q_dot_0 + k_scaled*t;
                qdd_des = k_scaled;
            else
                % compute trajectory at peak
                q_peak = q_0 + q_dot_0*t_plan + (1/2)*k_scaled*t_plan^2;
                qd_peak = q_dot_0 + k_scaled*t_plan;

                % compute second half of trajectory
                q_des = q_peak + qd_peak*(t-t_plan) + (1/2)*((0 - qd_peak)/(t_stop - t_plan))*(t-t_plan)^2;
                qd_des = qd_peak + ((0 - qd_peak)/(t_stop - t_plan))*(t-t_plan);
                qdd_des = (0 - qd_peak)/(t_stop - t_plan);
            end
        else
            % bring the trajectory to a stop in t_plan seconds
            % trajectory peaks at q_0
            q_peak = q_0;
            qd_peak = q_dot_0;
            
            if t <= t_plan && ~all(q_dot_0 == 0) % we're braking!
                q_des = q_peak + qd_peak*t + (1/2)*((0 - qd_peak)/t_plan)*t^2;
                qd_des = qd_peak + ((0 - qd_peak)/t_plan)*t;
                qdd_des = (0 - qd_peak)/t_plan;
            else % we should already be stopped, maintain that.
                q_des = q_peak;
                qd_des = zeros(size(q_dot_0));
                qdd_des = zeros(size(q_0));
            end
        end
    case 'bernstein'
        % assuming K = [-1, 1] corresponds to final position for now!!
        n_q = length(q_0);
        if ~isnan(k)
            q1 = q_0 + P.jrs_info.c_k_bernstein + P.jrs_info.g_k_bernstein.*k;
            for j = 1:n_q
                beta{j} = match_deg5_bernstein_coefficients({q_0(j); q_dot_0(j); q_ddot_0(j); q1(j); 0; 0});
                alpha{j} = bernstein_to_poly(beta{j}, 5);
            end
            q_des = zeros(length(q_0), 1);
            qd_des = zeros(length(q_0), 1);
            qdd_des = zeros(length(q_0), 1);
            for j = 1:n_q
                for coeff_idx = 0:5
                    q_des(j) = q_des(j) + alpha{j}{coeff_idx+1}*t^coeff_idx;
                    if coeff_idx > 0
                        qd_des(j) = qd_des(j) + coeff_idx*alpha{j}{coeff_idx+1}*t^(coeff_idx-1);
                    end
                    if coeff_idx > 1
                        qdd_des(j) = qdd_des(j) + (coeff_idx)*(coeff_idx-1)*alpha{j}{coeff_idx+1}*t^(coeff_idx-2);
                    end
                end
            end
        else
            % bring the trajectory to a stop using previous trajectory...
            t_plan = P.t_plan;
            if t <= t_plan && norm(q_dot_0) > 1e-8 && norm(q_dot_0) > 1e-8
                % just plug into previous trajectory, but shift time forward by t_plan.
                [q_des, qd_des, qdd_des] = P.info.desired_trajectory{end - 1}(t + t_plan);
            else % we should already be stopped, maintain that.
                q_des = q_0;
                qd_des = zeros(n_q, 1);
                qdd_des = zeros(n_q, 1);
            end
        end
    otherwise
        error('trajectory type not recognized');
    end
end

function [xp,yp] = circle(d,x,y,ang_start,ang_end,step)
    ang=ang_start:step:ang_end; 
    xp=(d/2)*cos(ang);
    yp=(d/2)*sin(ang);
end

function make_animation( h,index,filename )
    drawnow
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if index == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.001);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.001);
    end
end

function [M, C, g] = calculate_dynamics(q, qd, params)
    M = rnea_mass(q, params);
    C = rnea_coriolis(q, qd, params);
    g = rnea_gravity(q, params);
end