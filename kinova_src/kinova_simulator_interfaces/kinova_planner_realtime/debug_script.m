% 0. provide proper input data in buffer/armour.in
% 1. run ./compile_debug_script.sh
% 2. run ./test

close all; clear; clc;

%% initialize robot
robot = importrobot('Kinova_Grasp_URDF.urdf');
robot.DataFormat = 'col';
robot.Gravity = [0 0 -9.81];
params = load_robot_params(robot, ...
                           'add_uncertainty_to', 'all', ...
                           'uncertain_mass_range', [0.97, 1.03]);

A = uarmtd_agent(robot, params,...
                 'animation_set_axes_flag', 0,... 
                 'animation_set_view_flag', 0,...
                 'use_CAD_flag', true,...
                 'add_measurement_noise_', false, ...
                 'measurement_noise_size_', 0);
transmision_inertia = [8.02999999999999936 11.99620246153036440 9.00254278617515169 11.58064393167063599 8.46650409179141228 8.85370693737424297 8.85873036646853151]; % matlab doesn't import these from urdf so hard code into class

link_poly_zonotopes = create_pz_bounding_boxes(robot);

%% initialize desired trajectories
% choose random initial conditions and make sure they are aligned with
% the first three rows in buffer/armour.in
% q0 = [0.9534000000 -1.4310000000 0.1330000000 0.6418000000 -0.9534000000 -0.9534000000 0.0637000000 ]';
% qd0 = [0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000]'; 
% qdd0 = [0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000]';
q0 = [0.0000000000 -1.5707963268 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 ]';
qd0 = [1.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 ]';
qdd0 = [1.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 ]';
qdes = [0.4000000000 -1.5707963268 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 ]';

% choose a random k_range and make sure they are aligned with k_range in
% Parameters.h
% k_range = [pi/24, pi/24, pi/24, pi/24, pi/24, pi/24, pi/24]';
k_range = [pi/72, pi/72, pi/72, pi/72, pi/72, pi/72, pi/72]';


% choose a random point to slice and make sure they are equal to variable
% factors defined in PZ_test.cpp

k = [0.5, 0.7, 0.7, 0.0, -0.8, -0.6, -0.7]';
% k = zeros(7,1);
% k = ones(7,1);
% k = -ones(7,1);
% k=[0.999901, -0.898396, 0.973896, -0.454472, 0.999083, 0.940592, 0.974109]'

q1 = q0 + k .* k_range;
qd1 = zeros(7,1);
qdd1 = zeros(7,1);

beta = match_deg5_bernstein_coefficients({q0, qd0, qdd0, q1, qd1, qdd1});

% for tid = 1:128
duration = 2;
tid = 100;
tspan = linspace(0, duration, tid + 1);

%% read CUDA output

link_reachset_center = readmatrix('buffer/armour_joint_position_center.out', 'FileType', 'text');
link_reachset_generators = readmatrix('buffer/armour_joint_position_radius.out', 'FileType', 'text');

torque_reachset_center = readmatrix('buffer/armour_constraints.out', 'FileType', 'text');
torque_reachset_radius = readmatrix('buffer/armour_control_input_radius.out', 'FileType', 'text');

force_reachset_values = readmatrix('buffer/armour_wrench_values.out', 'FileType', 'text');
force_constraint_values = readmatrix('buffer/armour_force_constraint_radius.out', 'FileType', 'text');

des_traj_slice = readmatrix('buffer/armour_desired_sliced.out','FileType','text');

%% Processing CUDA output

% separate the desired trajectories
des_vel_center = des_traj_slice(:,1:2:14);
des_vel_radius = des_traj_slice(:,2:2:14);
des_aux_vel_center = des_traj_slice(:,15:2:28);
des_aux_vel_radius = des_traj_slice(:,16:2:28);
des_accel_center = des_traj_slice(:,29:2:end);
des_accel_radius = des_traj_slice(:,30:2:end);

% separate the force arrays
f_rs_c = force_reachset_values(:,1:3);
n_rs_c = force_reachset_values(:,4:6);
f_rs_r = force_reachset_values(:,7:9);
n_rs_r = force_reachset_values(:,10:12);
sep_ub_cuda = force_constraint_values(1:100,1);
slip_ub_cuda = force_constraint_values(101:200,1);
tip_ub_cuda = force_constraint_values(201:300,1);
sep_lb_cuda = force_constraint_values(1:100,2);
slip_lb_cuda = force_constraint_values(101:200,2);
tip_lb_cuda = force_constraint_values(201:300,2);

%% Verification

%% Plotting Link Reach Sets

figure; view(3); axis equal; hold on; axis on;

% choose a random time inside this time interval
t_lb = tspan(tid);
t_ub = tspan(tid + 1);
t = (t_ub - t_lb) * rand + t_lb;

q = get_desired_traj(beta, t, duration);

figure(1)
% plot robot
A.plot_at_time(q);

% ! Ask Bohao about the number 10 vs 7 in the for loop below
% plot link reachsets
numBodies = 8;
% for j = 1:robot.NumBodies
%     c = link_reachset_center((tid-1)*numBodies+j, :)';
%     g = link_reachset_generators( ((tid-1)*numBodies+j-1)*3+1 : ((tid-0)*numBodies+j)*3, :);
%     Z = zonotope(c, g);
%     Z_v = vertices(Z)';
%     trisurf(convhulln(Z_v),Z_v(:,1),Z_v(:,2),Z_v(:,3),'FaceColor',[0,0,1],'FaceAlpha',0.1,'EdgeColor',[0,0,1],'EdgeAlpha',0.3);
% end
% lighting flat
% end

%% Calculating Nominal Values

figure; hold on;

us = zeros(7,tid);
fs = cell(1,tid);
ns = cell(1,tid);
ts = zeros(1,tid);
for i = 1:tid
    % choose a random time inside this time interval
    t_lb = tspan(i);
    t_ub = tspan(i + 1);
    ts(i) = (t_ub - t_lb) * rand + t_lb;

    [q, qd, qdd] = get_desired_traj(beta, ts(i), duration);
    
    q_des_matlab(:,i) = q;
    qd_des_matlab(:,i) = qd;
    qdd_des_matlab(:,i) = qdd;

    [us(:,i), fs{i}, ns{i}] = rnea(q, qd, qd, qdd, true, params.nominal); % + transmision_inertia' .* qdd;
end

%% Plotting Torque Reach Sets

% u_lb = torque_reachset_center - torque_reachset_radius;
% u_ub = torque_reachset_center + torque_reachset_radius;
% 
% figure(2)
% % there is a better way to do this
% for i = 1:7
%     subplot(3,3,i);
%     hold on;
%     plot(ts, us(i,:), 'r');
%     plot(ts, u_lb(:,i), 'b');
%     plot(ts, u_ub(:,i), 'b');
%     title(['link ', num2str(i)]);
%     xlabel('time (sec)');
%     ylabel('torque (N*m)');
% end
% sgtitle('sliced torque reachable set');

%% Plotting Force and Moment Reach Sets

% extract the nominal values from the cells
for j = 1:tid
    f_nom(j,:) = fs{j}(:,10)';
    n_nom(j,:) = ns{j}(:,10)';
end

% get the upper and lower bounds of the force reach sets
f_ub = f_rs_c + f_rs_r;
f_lb = f_rs_c - f_rs_r;
n_ub = n_rs_c + n_rs_r;
n_lb = n_rs_c - n_rs_r;

figure(3)
plot_label = {'X-axis','Y-axis','Z-axis'};
for i = 1:3
    subplot(3,2,i*2-1)
    hold on
    plot(ts,f_nom(:,i),'ok')
    plot(ts,f_ub(:,i),'-b')
    plot(ts,f_lb(:,i),'-b')
    title([plot_label(i),' Force'])
    xlabel('Time (sec)')
    ylabel('Force (Newton)')
end
for i = 1:3
    subplot(3,2,i*2)
    hold on
    plot(ts,n_nom(:,i),'ok')
    plot(ts,n_ub(:,i),'-b')
    plot(ts,n_lb(:,i),'-b')
    title([plot_label(i),' Moment'])
    xlabel('Time (sec)')
    ylabel('Moment (Newton*meter)')
end

%% Calculate the Constraints

u_s = 0.609382421;
surf_rad =  0.058/2;

for i = 1:tid
    % separation constraint
    con_mat(i,1) = -1*f_nom(i,3);
    % slip constraint
    con_mat(i,2) = f_nom(i,1)^2 + f_nom(i,2)^2 - u_s^2*f_nom(i,3)^2;
    % tip constraint
    ZMP_top = cross([0;0;1],n_nom(i,:));
    ZMP_bottom = dot([0;0;1],f_nom(i,:));
    con_mat(i,3) = ZMP_top(1)^2+ZMP_top(2)^2 - surf_rad^2*ZMP_bottom^2;
end

figure(4)
constraint_label = {'Separation Constraint','Slipping Constraint','Tipping Constraint'};
for i = 1:3
    subplot(3,1,i)
    hold on
    plot(ts,con_mat(:,i),'-k')
    plot(ts,force_constraint_values((1+(i-1)*100):(100+(i-1)*100),1),'b-')
    plot(ts,force_constraint_values((1+(i-1)*100):(100+(i-1)*100),2),'b-')
    title(constraint_label(i))
    xlabel('Time (sec)')
    ylabel('Constraint Value')
end

%% Plotting the Force Constraints

% figure(4)
% for i=1:3
% 
% end

% for i = 1:length(A.time)
% 
%     out = W.grasp_check(A,A.agent_info,P.info)
% 
% end

%% Plotting Desired Trajectory Comparison

figure(5)
clf(5)
title('desired comparison')
for i = 1:7
    subplot(7,1,i)
    hold on
    plot(ts,qd_des_matlab(i,:),'ok')
    plot(ts,des_vel_center(:,i),'--r')
    plot(ts,des_aux_vel_center(:,i),'--b')
    % plot(ts,des_vel_center+des_vel_radius,'--r')
    % plot(ts,des_vel_center-des_vel_radius,'--r')
end

figure(6)
for i = 1:7
    subplot(7,1,i)
    hold on
    plot(ts,qdd_des_matlab(i,:),'ok')
    plot(ts,des_accel_center(:,i),'--b')
    % plot(ts,des_accel_center+des_vel_radius,'--r')
    % plot(ts,des_accel_center-des_vel_radius,'--r')
end


%% helper functions
function [q, qd, qdd] = get_desired_traj(beta, t, duration)
    [B, dB, ddB] = Bezier_kernel_deg5(t/duration); %t/dur
    
    q = zeros(7,1);
    qd = zeros(7,1);
    qdd = zeros(7,1);
    for j = 1:6
        q = q + beta{j} * B(j);
        qd = qd + 1/duration*beta{j} * dB(j);
        qdd = qdd + (1/duration)^2*beta{j} * ddB(j); % coeff
    end
end
