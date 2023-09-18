%% Kinova 3D Friction Cone Plotting

clear; 
close all; 
clc;

%%% To Do
% plot the nominal only with better shading
% plot the 


%% Setup Robot

% urdf names
% robot_name = 'Kinova_Grasp_URDF.urdf';
robot_name = 'Kinova_Grasp_Cylinder_Edge.urdf';

% import robot
robot = importrobot(robot_name);
robot.DataFormat = 'col';
robot.Gravity = [0 0 -9.81];

% get params
add_uncertainty_to = 'all'; % choose 'all', 'link', or 'none'
links_with_uncertainty = {}; % if add_uncertainty_to = 'link', specify links here.
uncertain_mass_range = [0.97, 1.03];
params = load_robot_params(robot, ...
                           'add_uncertainty_to', add_uncertainty_to, ...
                           'links_with_uncertainty', links_with_uncertainty,...
                           'uncertain_mass_range', uncertain_mass_range);

% controller info
LLC_info.ultimate_bound = 0.00191; % 0.00191
LLC_info.Kr = 10;

% create link poly zonotopes
[link_poly_zonotopes, link_sizes, meshes] = create_link_poly_zonos(robot);
zono_order = 10;

% contact parameters
u_s = 0.5;
surf_rad = 0.05; % meter

%% Setup Plotting Flags

plot_idx = 0;

plot_force_trajectory = false;
plot_trajectory_1 = false;

save_plot = false;

%% Color Coding

% zonotope coloring
unsliced_color = 1/256*[90,200,243];
% slice_color = 1/256*[72,181,163]; % light green
slice_color = 1/256*[75,154,76]; % dark green
face_alpha = 0.3;
face_alpha_light = 0.03;

%% Plot Formatting

fontsize = 16;

%% compute trajectories

% initial conditions
% old figure
% q_0 = [0;-pi/4;0;-pi/2;0;pi/4;0];
% qd_0= [-pi/3;-pi/12;0;0;0;0;0];
% qdd_0 = [pi/3;-pi/6;0;0;0;0;-pi/24];
% new figure
q_0 = [pi/4;-pi/4;0;-pi/2;0;pi/4;0];
qd_0= [0;pi/30;0;0;0;0;0];
qdd_0 = [0;pi/30;0;0;0;0;0];

% plotting to verify configuration
% clf(101)
% figure(101)
% show(robot,q_0)

% create pz trajectories
joint_axes = [zeros(2, length(q_0)); ones(1, length(q_0))]; % Question: shouldn't this be pulled out of the urdf or params variable?
taylor_degree = 5;
traj_type = 'bernstein';
add_ultimate_bound = true;

% calculating JRS
[Q_des, Qd_des, Qdd_des, Q, Qd, Qd_a, Qdd_a, R_des, R_t_des, R, R_t, T, E_p, jrs_info] = create_jrs_online_modified(q_0,qd_0,qdd_0, joint_axes, taylor_degree, traj_type, add_ultimate_bound, LLC_info);

%% Trajectory Parameter

% trajectory parameter
% old figure
% kvec = [0.6; -0.8; 0.5; -0.2; -0.4; 0.35; 0.34];
% new figure
kvec = [0.001; 0.5; 0.001; 0.5; 0.001; 0.5; 0.001];
% kvec = [-1;-1;-1;-1;-1;-1;-1];
% kvec = [1;1;1;1;1;1;1];

%% Desired Trajectories

% planning info
t_traj = 0:jrs_info.dt:jrs_info.t_f;
t_traj(end) = [];

id_tmp = jrs_info.id;
id_names_tmp = jrs_info.id_names;
jrs_info.id = {};
jrs_info.id_names = {};
for i=1:length(kvec)
    jrs_info.id{i,1} = id_tmp(i,:);
    jrs_info.id_names{i,1} = id_names_tmp(i,:);
end

P.jrs_info = jrs_info;
P.t_start = 0;
P.t_stop = 1;
P.t_plan = 0.5;
P.traj_type = traj_type;

% compute desired trajectory
t_steps_dt = 0:jrs_info.dt:1;
t_steps = 0:jrs_info.dt:1;
n_steps = length(t_steps_dt);
q_des = zeros(length(kvec), length(t_steps));
q_des_dt = zeros(length(kvec), n_steps);
qd_des = zeros(length(kvec), n_steps);
qdd_des = zeros(length(kvec), n_steps);

for i = 1:length(t_steps)
    [q_des(:,i), qd_des(:,i), qdd_des(:,i)] = desired_trajectory(P, q_0, qd_0, qdd_0, t_steps(i), kvec);
end

for i = 1:jrs_info.n_t
    [q_des_dt(:,i), ~, ~] = desired_trajectory(P, q_0, qd_0, qdd_0, t_steps_dt(i), kvec);
end

% compute trajectory bounds
q_max = q_des + LLC_info.ultimate_bound / LLC_info.Kr;
q_min = q_des - LLC_info.ultimate_bound / LLC_info.Kr;

% ids for slicing
id_slice = 1:length(kvec);

% compute error zonotopes
Q_e = cell(jrs_info.n_t, 1);
R_e = cell(jrs_info.n_t, 1);
R_t_e = cell(jrs_info.n_t, 1);

for i = 1:jrs_info.n_t
	for j = 1:jrs_info.n_q
        Q_e{i}{j, 1} = q_des_dt(j,i) + E_p{j};
		[R_e{i}{j, 1}, R_t_e{i}{j, 1}] = get_pz_rotations_from_q(Q_e{i}{j, 1}, joint_axes(:, j), taylor_degree);
	end
end

%% Calling PZRNEA

for i = 1:jrs_info.n_t
    [tau_temp, f_temp, n_temp] = poly_zonotope_rnea(R{i}, R_t{i}, Qd{i}, Qd_a{i}, Qdd_a{i}, true, params.pz_interval);
%     tau_int{i} = tau_temp{10,1};
    f_int{i} = f_temp{10,1};
    n_int{i} = n_temp{10,1};
end

%% Calling RNEA for Nominal Wrench Trajectories

for i = 1:length(t_steps)
    [tau_temp, f_temp, n_temp] = rnea(q_des(:,i), qd_des(:,i), qd_des(:,i), qdd_des(:,i), true, params.nominal);
%     tau_int{i} = tau_temp{10,1};
    f_nom(:,i) = f_temp(:,10);
    n_nom(:,i) = n_temp(:,10);
end

%% Calculating the PZ Force Constraints

for i = 1:floor(length(t_traj)/2)

    % separation constraint
%     sep_PZ{1,i} = -1*f_int{1,i};
%     
    temp = f_int{1,i};
    f_int_x = polyZonotope_ROAHM(temp.c(1),temp.G(1,:),temp.Grest(1,:),temp.expMat,temp.id);
    f_int_y = polyZonotope_ROAHM(temp.c(2),temp.G(2,:),temp.Grest(2,:),temp.expMat,temp.id);
    f_int_z = polyZonotope_ROAHM(temp.c(3),temp.G(3,:),temp.Grest(3,:),temp.expMat,temp.id);
    
    % slipping constraint
    slip_PZ{1,i} = f_int_x.*f_int_x+f_int_y.*f_int_y - u_s^2*f_int_z.*f_int_z;

    % calculating the PZ form of the ZMP constraint
%     ZMP_PZ_top = cross([0;0;1],n_int{1,i});
%     ZMP_PZ_bottom = f_int{1,i}*[0,0,1];
%     ZMP_PZ_bottom_int = interval(ZMP_PZ_bottom);
%     ZMP_PZ_bottom_inf = ZMP_PZ_bottom_int.inf;
%     ZMP_PZ_bottom_sup = ZMP_PZ_bottom_int.sup;
%     % Note that this is not the entire overapproximated PZ because the
%     % bottom has the .sup as well
%     ZMP_PZ_inf = interval(ZMP_PZ_top) / ZMP_PZ_bottom_inf;
%     ZMP_PZ_sup = interval(ZMP_PZ_top) / ZMP_PZ_bottom_sup;
%     ZMP_PZ{i} = convHull(ZMP_PZ_inf,ZMP_PZ_sup);
% 
%     % calculating sliced version of the ZMP constraint overapproximation
%     % slice the ZMP_PZ_TOP and then do the division and take convHull?
%     % can also slice the bottom first too
%     ZMP_PZ_bottom_sliced = getSubset(ZMP_PZ_bottom,ZMP_PZ_bottom.id,kvec(ZMP_PZ_bottom.id));
%     ZMP_PZ_bottom_sliced_int = interval(ZMP_PZ_bottom_sliced);
%     ZMP_PZ_bottom_sliced_inf = ZMP_PZ_bottom_sliced_int.inf;
%     ZMP_PZ_bottom_sliced_sup = ZMP_PZ_bottom_sliced_int.sup;
% 
%     ZMP_PZ_top_sliced = getSubset(ZMP_PZ_top,ZMP_PZ_top.id,kvec(ZMP_PZ_top.id));
%     ZMP_PZ_sliced_inf = interval(ZMP_PZ_top_sliced) / ZMP_PZ_bottom_sliced_inf;
%     ZMP_PZ_sliced_sup = interval(ZMP_PZ_top_sliced) / ZMP_PZ_bottom_sliced_sup;
%     ZMP_PZ_sliced{i} = convHull(ZMP_PZ_sliced_inf,ZMP_PZ_sliced_sup);

end

%% Calculating the Nominal Constraints

% separation constraint
% sep_nom = -1*f_nom(3,:);

% slipping constraint
slip_nom = sqrt(f_nom(1,:).^2+f_nom(2,:).^2) - u_s.*abs(f_nom(3,:));
slip2_nom = f_nom(1,:).^2+f_nom(2,:).^2 - u_s^2.*f_nom(3,:).^2; % NOTE: probably want to use this one

% for i = 1:floor(length(t_steps)/2)
%     % tipping constraint the normal way
%     ZMP_top = cross([0;0;1],n_nom(:,i)); % normal vector should come first
%     ZMP_bottom = dot([0;0;1],f_nom(:,i));
%     ZMP(:,i) = ZMP_top/ZMP_bottom;
%     ZMP_rad(i) = sqrt(ZMP(1,i)^2+ZMP(2,i)^2);
%     tip(i) = ZMP_rad(i) - surf_rad;
%     % tipping constraint the PZ way
%     ZMP_top2 = cross([0;0;1],n_nom(:,i));
%     ZMP_bottom2 = dot([0;0;1],f_nom(:,i));
%     tip2(i) = ZMP_top(1)^2 + ZMP_top(2)^2 - ZMP_bottom^2*(surf_rad)^2;
% end

%% Plotting Wrench Trajectory

if plot_force_trajectory

    % choose which force to plot (1=x-axis,
    force = 3; 

    % plotting the unsliced force trajectory
    plot_idx = plot_idx + 1;
    figure(plot_idx); clf; hold on;
    title('Contact Joint z-Axis Force')

    for i = 1:floor(length(t_traj)/2)
        % plot the polynomial overapproximation
        % calculate the inf/sup
        if f_int{1,i}.G
            poly_inf = f_int{1,i}.c(force) - sum(abs(f_int{1,i}.G(force,:))) - sum(abs(f_int{1,i}.Grest(force,:)));
            poly_sup = f_int{1,i}.c(force) + sum(abs(f_int{1,i}.G(force,:))) + sum(abs(f_int{1,i}.Grest(force,:)));
        else % the magnitude of Grest means that only those were tracked
            poly_inf = f_int{1,i}.c(force) - sum(abs(f_int{1,i}.Grest(force,:)));
            poly_sup = f_int{1,i}.c(force) + sum(abs(f_int{1,i}.Grest(force,:)));
        end
        p1 = patch([t_traj(i)+jrs_info.dt; t_traj(i)+jrs_info.dt; t_traj(i); t_traj(i)], [poly_sup; poly_inf; poly_inf; poly_sup],'b');
%         p1.EdgeColor = pz_err_color;
        p1.LineWidth = 0.1;
        p1.FaceColor = unsliced_color;
        p1.FaceAlpha = face_alpha;

    end

    % plot the nominal values
    plot(t_steps, f_nom(force,:),'-k')

    % formatting for plot
    ylim([0 2.25])

    % plotting the sliced force trajectory
%     plot_idx = plot_idx + 1;
%     figure(plot_idx); clf; hold on;
%     title('Force Plot: Last Joint')

    for i = 1:floor(length(t_traj)/2)
        f_sliced{1,i} = getSubset(f_int{1,i}, f_int{1,i}.id, kvec(f_int{1,i}.id));
        if f_sliced{1,i}.G
            poly_inf = f_sliced{1,i}.c(force) - sum(abs(f_sliced{1,i}.G(force,:))) - sum(abs(f_sliced{1,i}.Grest(force,:)));
            poly_sup = f_sliced{1,i}.c(force) + sum(abs(f_sliced{1,i}.G(force,:))) + sum(abs(f_sliced{1,i}.Grest(force,:)));
        else
            poly_inf = f_sliced{1,i}.c(force) - sum(abs(f_sliced{1,i}.Grest(force,:)));
            poly_sup = f_sliced{1,i}.c(force) + sum(abs(f_sliced{1,i}.Grest(force,:)));
        end
        p1 = patch([t_traj(i)+jrs_info.dt; t_traj(i)+jrs_info.dt; t_traj(i); t_traj(i)], [poly_sup; poly_inf; poly_inf; poly_sup],'b');
%         p1.EdgeColor = pz_err_color;
        p1.LineWidth = 0.1;
        p1.FaceColor = slice_color;
        p1.FaceAlpha = face_alpha;
    end
    % plot the nominal values
    plot(t_steps, f_nom(force,:),'-k')

    % formatting for plot
    ylim([0 2.25])
    xlabel('Time (s)')
    ylabel('Force (N)')

    set(gca,'FontSize',fontsize)

end

%% Plotting Nominal Friction Cone Only

plot_idx = plot_idx + 1;
figure(plot_idx); clf; % hold on;
subplot(1,2,2)
hold on

theta = linspace(0,2*pi,50);

% plot the friction cone boundary
max_fz = max(f_nom(3,:));
min_fz = min(f_nom(3,:));

plotConcentricCircles(plot_idx,min_fz*u_s,max_fz*u_s,[1,0,0,0.25],'r')

% plotting the nominal values
num = length(t_steps); % floor(length(t_steps)/2);
plot(f_nom(1,1:num),f_nom(2,1:num),'-k','LineWidth',3)

%%% formatting
xlabel('x-axis Tangential Force (N)')
ylabel('y-axis Tangential Force (N)')
xticks([-1 -0.5 0 0.5 1])
yticks([-1 -0.5 0 0.5 1])
axis('square')
% grid on
box on
set(gca,'FontSize',fontsize)
set(gcf, 'Color', 'w');

% Plotting Nominal 3D Friction Cone Only

% plot_idx = plot_idx + 1;
% figure(plot_idx); clf; hold on;
subplot(1,2,1)
hold on

%%% calculating continuous force trajectory
% time
t_start = t_traj(1);
% t_start = 10/40;
% t_end = 1/40; % for single time interval
% t_end = 11/40; % for full planning iteration
t_end = t_traj(floor(length(t_traj)/2));
t_cont = linspace(t_start,t_end);
% desired trajectory
for i = 1:length(t_cont)
    [q_cont_des(:,i), qd_cont_des(:,i), qdd_cont_des(:,i)] = desired_trajectory(P, q_0, qd_0, qdd_0, t_cont(i), kvec);
    % rnea
    [u_temp, f_temp, n_temp] = rnea(q_cont_des(:,i), qd_cont_des(:,i), qd_cont_des(:,i), qdd_cont_des(:,i), true, params.nominal);
    f_cont(:,i) = f_temp(:,10);
    n_cont(:,i) = n_temp(:,10);
end

%%% plot the friction cone slice boundaries
theta_friction = linspace(0,2*pi,100);
% plot rings for all time steps
% for i = 1:length(t_steps)
%     r_friction = f_nom(3,i)*u_s;
%     plot3(r_friction*cos(theta_friction),r_friction*sin(theta_friction),f_nom(3,i)*ones(1,length(theta_friction)),'-r')
% end
% plot rings for max and minimum normal force
min_r_friction = min(f_nom(3,1:floor(length(t_traj)/2)))*u_s;
max_r_friction = max(f_nom(3,1:floor(length(t_traj)/2)))*u_s;
plot3(min_r_friction*cos(theta_friction),min_r_friction*sin(theta_friction),min(f_nom(3,1:floor(length(t_traj)/2)))*ones(1,length(theta_friction)),'-r')
plot3(max_r_friction*cos(theta_friction),max_r_friction*sin(theta_friction),max(f_nom(3,1:floor(length(t_traj)/2)))*ones(1,length(theta_friction)),'-r')

%%% plot the shaded friction cone
% active portion
r = linspace(min_r_friction,max_r_friction,10);
theta = linspace(0,2*pi,50);
[RR,Theta] = meshgrid(r,theta);
X = RR.*cos(Theta);
Y = RR.*sin(Theta);
Z = 1./u_s.*RR;
h1 = surf(X,Y,Z,'EdgeColor','none','FaceColor','r','FaceAlpha',0.25);
% below active portion
r = linspace(0,min_r_friction,10);
theta = linspace(0,2*pi,50);
[RR,Theta] = meshgrid(r,theta);
X = RR.*cos(Theta);
Y = RR.*sin(Theta);
Z = 1./u_s.*RR;
h1 = surf(X,Y,Z,'EdgeColor','none','FaceColor','r','FaceAlpha',0.05);
% above active portion
r = linspace(max_r_friction,1,10);
theta = linspace(0,2*pi,50);
[RR,Theta] = meshgrid(r,theta);
X = RR.*cos(Theta);
Y = RR.*sin(Theta);
Z = 1./u_s.*RR;
h1 = surf(X,Y,Z,'EdgeColor','none','FaceColor','r','FaceAlpha',0.05);

%%% plot the nominal trajectory
plot3(f_cont(1,:),f_cont(2,:),f_cont(3,:),'-k', 'LineWidth', 3)

%%% formatting
xlabel('x-axis Tangential Force (N)')
ylabel('y-axis Tangential Force (N)')
zlabel('z-axis Normal Force (N)')
xticks([-1 -0.5 0 0.5 1])
yticks([-1 -0.5 0 0.5 1])
axis('square')
% grid on
box on
set(gca,'FontSize',fontsize)
set(gcf, 'Color', 'w');
view(-45,20)

%% Plotting Single Time Interval Friction Cone

plot_idx = plot_idx + 1;
figure(plot_idx); clf; % hold on;
% title('Friction Cone Unsliced Plot')
% title('Friction Cone')
subplot(1,2,2)
hold on

theta = linspace(0,2*pi,50);

% plot the friction cone boundary
max_fz = max(f_nom(3,:));
min_fz = min(f_nom(3,:));

plotConcentricCircles(plot_idx,min_fz*u_s,max_fz*u_s,[1,0,0,0.25],'r')

% plot the unsliced force overapproximation
for i = 1:1 % floor(length(t_traj)/2)
    f_int_zono_2D = zonotope(f_int{i}); % not needed?
    fp1 = plot(f_int_zono_2D,[1,2],'Filled',true,'FaceColor',unsliced_color);
    fp1.LineWidth = 0.1;
    fp1.FaceColor = unsliced_color;
    fp1.EdgeColor = unsliced_color;
    fp1.FaceAlpha = face_alpha_light;
end

% for separate 2D friction cone plots
% plot_idx = plot_idx + 1;
% figure(plot_idx); clf; hold on;
% title('Friction Cone Sliced Plot')

% plot sliced force overapproximation
for i = 1:1 % floor(length(t_traj)/2)
    % slice
    f_sliced{i} = getSubset(f_int{i},f_int{i}.id,kvec(f_int{i}.id));
    % plot
    fp2 = plot(interval(f_sliced{i}),[1,2],'Filled',true);
    fp2.LineWidth = 0.1;
    fp2.FaceColor = slice_color;
    fp2.EdgeColor = slice_color;
    fp2.FaceAlpha = face_alpha;
end

% plotting the nominal values
% num = floor(length(t_steps)/40); % length(t_steps); % 
plot(f_nom(1,1:2),f_nom(2,1:2),'-k','LineWidth',3)

xlabel('x-axis Force (N)')
ylabel('y-axis Force (N)')
set(gca,'FontSize',fontsize)
set(gcf, 'Color', 'w');
axis square
% axis equal
xticks([-1 -0.5 0 0.5 1])
yticks([-1 -0.5 0 0.5 1])
box on

% Plotting Friction Cone and Force PZ in 3D

% plot_idx = plot_idx + 1;
% figure(plot_idx); clf; hold on;
% title('Friction Cone 3D Plot')
subplot(1,2,1)
hold on

%%% calculating continuous force trajectory
% time
t_start = t_traj(1);
% t_start = 10/40;
% t_end = 1/40; % for single time interval
% t_end = 11/40; % for full planning iteration
t_end = t_traj(2);
t_cont = linspace(t_start,t_end);
% desired trajectory
for i = 1:length(t_cont)
    [q_cont_des(:,i), qd_cont_des(:,i), qdd_cont_des(:,i)] = desired_trajectory(P, q_0, qd_0, qdd_0, t_cont(i), kvec);
    % rnea
    [u_temp, f_temp, n_temp] = rnea(q_cont_des(:,i), qd_cont_des(:,i), qd_cont_des(:,i), qdd_cont_des(:,i), true, params.nominal);
    f_cont(:,i) = f_temp(:,10);
    n_cont(:,i) = n_temp(:,10);
end

%%% plot the shaded friction cone
% active portion
r = linspace(min_r_friction,max_r_friction,10);
theta = linspace(0,2*pi,50);
[RR,Theta] = meshgrid(r,theta);
X = RR.*cos(Theta);
Y = RR.*sin(Theta);
Z = 1./u_s.*RR;
h1 = surf(X,Y,Z,'EdgeColor','none','FaceColor','r','FaceAlpha',0.25);
% below active portion
r = linspace(0,min_r_friction,10);
theta = linspace(0,2*pi,50);
[RR,Theta] = meshgrid(r,theta);
X = RR.*cos(Theta);
Y = RR.*sin(Theta);
Z = 1./u_s.*RR;
h1 = surf(X,Y,Z,'EdgeColor','none','FaceColor','r','FaceAlpha',0.05);
% above active portion
r = linspace(max_r_friction,1,10);
theta = linspace(0,2*pi,50);
[RR,Theta] = meshgrid(r,theta);
X = RR.*cos(Theta);
Y = RR.*sin(Theta);
Z = 1./u_s.*RR;
h1 = surf(X,Y,Z,'EdgeColor','none','FaceColor','r','FaceAlpha',0.05);

%%% plot the friction cone slices
theta_friction = linspace(0,2*pi,100);
% plot rings for all time steps
% for i = 1:length(t_steps)
%     r_friction = f_nom(3,i)*u_s;
%     plot3(r_friction*cos(theta_friction),r_friction*sin(theta_friction),f_nom(3,i)*ones(1,length(theta_friction)),'-r')
% end
% plot rings for max and minimum normal force
min_r_friction = min(f_nom(3,1:floor(length(t_traj)/2)))*u_s;
max_r_friction = max(f_nom(3,1:floor(length(t_traj)/2)))*u_s;
plot3(min_r_friction*cos(theta_friction),min_r_friction*sin(theta_friction),min(f_nom(3,1:floor(length(t_traj)/2)))*ones(1,length(theta_friction)),'-r')
plot3(max_r_friction*cos(theta_friction),max_r_friction*sin(theta_friction),max(f_nom(3,1:floor(length(t_traj)/2)))*ones(1,length(theta_friction)),'-r')

%%% plot the nominal trajectory
plot3(f_cont(1,:),f_cont(2,:),f_cont(3,:),'-k', 'LineWidth', 3)

%%% plot the overapproximations
for i = 1:1:1 % floor(length(t_traj)/2)

    % testing reduction operations
%     plot(f_int{i}, [1,2,3]) % ,'Filled',true)
%     f_int_zono = zonotope(f_int{i});
%     plot(f_int_zono,[1,2,3],'EdgeColor','r')
%     f_int_zono_reduce = reduce(f_int_zono,'girard',1);
%     plot(f_int_zono_reduce,[1,2,3],'EdgeColor','c')

    %%% plot zonotope overapproximation form
    % plot the unsliced overapproximation
    f_int_zono = zonotope(f_int{i});
    f_int_reduce = reduce(f_int_zono,'girard',1);
    V_unsliced = vertices(f_int_reduce)';
    [V_convhull_unsliced, V_slc_unsliced] = convhull(V_unsliced(:,1),V_unsliced(:,2),V_unsliced(:,3));
    trisurf(V_convhull_unsliced,V_unsliced(:,1),V_unsliced(:,2),V_unsliced(:,3),'FaceColor',unsliced_color,'FaceAlpha',0.3,'EdgeAlpha',0.1)

    % plot the sliced overapproximation
    f_sliced_2 = getSubset(f_int{i},f_int{i}.id,kvec(f_int{i}.id));
%     f_int_convHull = convHull(f_sliced_2);
    f_int_zono = zonotope(f_sliced_2);
    f_int_reduce = reduce(f_int_zono,'girard',1);
    V_sliced = vertices(f_int_reduce)';
    [V_convhull_sliced, V_slc_sliced] = convhull(V_sliced(:,1),V_sliced(:,2),V_sliced(:,3));
    trisurf(V_convhull_sliced,V_sliced(:,1),V_sliced(:,2),V_sliced(:,3),'FaceColor',slice_color,'FaceAlpha',0.3,'EdgeAlpha',0.1)

end

%%% format plot
xlabel('x-axis Tangential Force (N)')
ylabel('y-axis Tangential Force (N)')
zlabel('z-axis Normal Force (N)')
axis('square')
% axis equal
set(gca,'FontSize',fontsize)
set(gcf, 'Color', 'w');
xticks([-1 -0.5 0 0.5 1])
yticks([-1 -0.5 0 0.5 1])
% grid on
box on
view(-45,20)

%% Plotting Full Planning Iteration Friction Cone

plot_idx = plot_idx + 1;
figure(plot_idx); clf; % hold on;
% title('Friction Cone Unsliced Plot')
% title('Friction Cone')
subplot(1,2,2)
hold on

theta = linspace(0,2*pi,50);

% plot the friction cone boundary
max_fz = max(f_nom(3,:));
min_fz = min(f_nom(3,:));

plotConcentricCircles(plot_idx,min_fz*u_s,max_fz*u_s,[1,0,0,0.25],'r')

% plot the unsliced force overapproximation
for i = 1:floor(length(t_traj)/2)
    f_int_zono_2D = zonotope(f_int{i}); % not needed?
    fp1 = plot(f_int_zono_2D,[1,2],'Filled',true,'FaceColor',unsliced_color);
    fp1.LineWidth = 0.1;
    fp1.FaceColor = unsliced_color;
    fp1.EdgeColor = unsliced_color;
    fp1.FaceAlpha = face_alpha_light;
end

% for separate 2D friction cone plots
% plot_idx = plot_idx + 1;
% figure(plot_idx); clf; hold on;
% title('Friction Cone Sliced Plot')

% plot sliced force overapproximation
for i = 1:floor(length(t_traj)/2)
    % slice
    f_sliced{i} = getSubset(f_int{i},f_int{i}.id,kvec(f_int{i}.id));
    % plot
    fp2 = plot(interval(f_sliced{i}),[1,2],'Filled',true);
    fp2.LineWidth = 0.1;
    fp2.FaceColor = slice_color;
    fp2.EdgeColor = slice_color;
    fp2.FaceAlpha = face_alpha;
end

% %%% calculating continuous force trajectory
% % time
% % t_start = 0;
% t_start = 10/40;
% % t_end = 1/40; % for single time interval
% t_end = 11/40; % for full planning iteration
% t_cont = linspace(t_start,t_end);
% % desired trajectory
% for i = 1:length(t_cont)
%     [q_cont_des(:,i), qd_cont_des(:,i), qdd_cont_des(:,i)] = desired_trajectory(P, q_0, qd_0, qdd_0, t_cont(i), kvec);
%     % rnea
%     [u_temp, f_temp, n_temp] = rnea(q_cont_des(:,i), qd_cont_des(:,i), qd_cont_des(:,i), qdd_cont_des(:,i), true, params.nominal);
%     f_cont(:,i) = f_temp(:,10);
%     n_cont(:,i) = n_temp(:,10);
% end

% plotting the nominal values
num = floor(length(t_steps)/2); % length(t_steps); % 
plot(f_nom(1,1:num),f_nom(2,1:num),'-k','LineWidth',3)

%%% format plot
xlabel('x-axis Tangential Force (N)')
ylabel('y-axis Tangential Force (N)')
axis('square')
% axis equal
set(gca,'FontSize',fontsize)
set(gcf, 'Color', 'w');
xticks([-1 -0.5 0 0.5 1])
yticks([-1 -0.5 0 0.5 1])
% grid on
box on

% Plotting Friction Cone and Force PZ in 3D

% plot_idx = plot_idx + 1;
% figure(plot_idx); clf; hold on;
% title('Friction Cone 3D Plot')
subplot(1,2,1)
hold on

%%% calculating continuous force trajectory
% time
t_start = t_traj(1);
% t_start = 10/40;
% t_end = 1/40; % for single time interval
% t_end = 11/40; % for full planning iteration
t_end = t_traj(floor(length(t_traj)/2));
t_cont = linspace(t_start,t_end);
% desired trajectory
for i = 1:length(t_cont)
    [q_cont_des(:,i), qd_cont_des(:,i), qdd_cont_des(:,i)] = desired_trajectory(P, q_0, qd_0, qdd_0, t_cont(i), kvec);
    % rnea
    [u_temp, f_temp, n_temp] = rnea(q_cont_des(:,i), qd_cont_des(:,i), qd_cont_des(:,i), qdd_cont_des(:,i), true, params.nominal);
    f_cont(:,i) = f_temp(:,10);
    n_cont(:,i) = n_temp(:,10);
end

%%% plot the shaded friction cone
% active portion
r = linspace(min_r_friction,max_r_friction,10);
theta = linspace(0,2*pi,50);
[RR,Theta] = meshgrid(r,theta);
X = RR.*cos(Theta);
Y = RR.*sin(Theta);
Z = 1./u_s.*RR;
h1 = surf(X,Y,Z,'EdgeColor','none','FaceColor','r','FaceAlpha',0.25);
% below active portion
r = linspace(0,min_r_friction,10);
theta = linspace(0,2*pi,50);
[RR,Theta] = meshgrid(r,theta);
X = RR.*cos(Theta);
Y = RR.*sin(Theta);
Z = 1./u_s.*RR;
h1 = surf(X,Y,Z,'EdgeColor','none','FaceColor','r','FaceAlpha',0.05);
% above active portion
r = linspace(max_r_friction,1,10);
theta = linspace(0,2*pi,50);
[RR,Theta] = meshgrid(r,theta);
X = RR.*cos(Theta);
Y = RR.*sin(Theta);
Z = 1./u_s.*RR;
h1 = surf(X,Y,Z,'EdgeColor','none','FaceColor','r','FaceAlpha',0.05);

%%% plot the friction cone slices
theta_friction = linspace(0,2*pi,100);
% plot rings for all time steps
% for i = 1:length(t_steps)
%     r_friction = f_nom(3,i)*u_s;
%     plot3(r_friction*cos(theta_friction),r_friction*sin(theta_friction),f_nom(3,i)*ones(1,length(theta_friction)),'-r')
% end
% plot rings for max and minimum normal force
min_r_friction = min(f_nom(3,1:floor(length(t_traj)/2)))*u_s;
max_r_friction = max(f_nom(3,1:floor(length(t_traj)/2)))*u_s;
plot3(min_r_friction*cos(theta_friction),min_r_friction*sin(theta_friction),min(f_nom(3,1:floor(length(t_traj)/2)))*ones(1,length(theta_friction)),'-r')
plot3(max_r_friction*cos(theta_friction),max_r_friction*sin(theta_friction),max(f_nom(3,1:floor(length(t_traj)/2)))*ones(1,length(theta_friction)),'-r')

%%% plot the nominal trajectory
plot3(f_cont(1,:),f_cont(2,:),f_cont(3,:),'-k', 'LineWidth', 3)

%%% plot the overapproximations
for i = 1:1:floor(length(t_traj)/2)

    % testing reduction operations
%     plot(f_int{i}, [1,2,3]) % ,'Filled',true)
%     f_int_zono = zonotope(f_int{i});
%     plot(f_int_zono,[1,2,3],'EdgeColor','r')
%     f_int_zono_reduce = reduce(f_int_zono,'girard',1);
%     plot(f_int_zono_reduce,[1,2,3],'EdgeColor','c')

    %%% plot zonotope overapproximation form
    % plot the unsliced overapproximation
    f_int_zono = zonotope(f_int{i});
    f_int_reduce = reduce(f_int_zono,'girard',1);
    V_unsliced = vertices(f_int_reduce)';
    [V_convhull_unsliced, V_slc_unsliced] = convhull(V_unsliced(:,1),V_unsliced(:,2),V_unsliced(:,3));
    trisurf(V_convhull_unsliced,V_unsliced(:,1),V_unsliced(:,2),V_unsliced(:,3),'FaceColor',unsliced_color,'FaceAlpha',0.03,'EdgeAlpha',0.05)

    % plot the sliced overapproximation
    f_sliced_2 = getSubset(f_int{i},f_int{i}.id,kvec(f_int{i}.id));
%     f_int_convHull = convHull(f_sliced_2);
    f_int_zono = zonotope(f_sliced_2);
    f_int_reduce = reduce(f_int_zono,'girard',1);
    V_sliced = vertices(f_int_reduce)';
    [V_convhull_sliced, V_slc_sliced] = convhull(V_sliced(:,1),V_sliced(:,2),V_sliced(:,3));
    trisurf(V_convhull_sliced,V_sliced(:,1),V_sliced(:,2),V_sliced(:,3),'FaceColor',slice_color,'FaceAlpha',0.3,'EdgeAlpha',0.1)

%     %%% plot PZ form
%     % plot the unsliced overapproximation
%     plot(f_int{i},[1,2,3],'Filled',true,'FaceColor','c','FaceAlpha',0.01)
% 
%     % plot the sliced overapproximation
%     f_sliced_2 = getSubset(f_int{i},f_int{i}.id,kvec(f_int{i}.id));
%     plot(f_sliced_2,[1,2,3],'Filled',true,'FaceColor','g','FaceAlpha',0.05)

end

%%% format plot
xlabel('x-axis Tangential Force (N)')
ylabel('y-axis Tangential Force (N)')
zlabel('z-axis Normal Force (N)')
axis('square')
% axis equal
set(gca,'FontSize',fontsize)
set(gcf, 'Color', 'w');
xticks([-1 -0.5 0 0.5 1])
yticks([-1 -0.5 0 0.5 1])
% grid on
box on
view(-45,20)

%% Plotting Animation of Friction Cone in 3D

% plot_idx = plot_idx + 1;
% fig1 = figure(plot_idx); 
% clf; hold on;
% title('Friction Cone 3D Plot')
% 
% for i = 1:1:length(t_traj)
% 
%     % plot the overapproximation
%     f_int_convHull = convHull(f_int{i});
%     f_int_zono = zonotope(f_int_convHull);
%     f_int_reduce = reduce(f_int_zono,'girard',1);
% %     fc1 = plot(f_int_reduce, [1,2,3]); %,[1,2,3],'Splits',1); %,'Filled',true);
% %     fc1.LineWidth = 0.1;
% %     fc1.FaceColor = slice_color;
% %     fc1.FaceAlpha = 0.3;
% %     fc1.EdgeAlpha = 0.3;
% 
%     % TODO: can plot the force PZ by getting the inf, sup for each and using fill3? or patch?
%     V = vertices(f_int_reduce)';
%     [V_convhull, V_slc] = convhull(V(:,1),V(:,2),V(:,3));
% %     trisurf(V_convhull,V(:,1),V(:,2),V(:,3),'FaceColor',unsliced_color,'FaceAlpha',0.03,'EdgeAlpha',0.1)
% 
%     % plot the sliced overapproximation
%     % TODO: plot the sliced overapproximation
%     f_sliced_2 = getSubset(f_int{i},f_int{i}.id,kvec(f_int{i}.id));
%     f_int_convHull = convHull(f_sliced_2);
%     f_int_zono = zonotope(f_int_convHull);
%     f_int_reduce = reduce(f_int_zono,'girard',1);
%     V = vertices(f_int_reduce)';
%     [V_convhull, V_slc] = convhull(V(:,1),V(:,2),V(:,3));
%     trisurf(V_convhull,V(:,1),V(:,2),V(:,3),'FaceColor',slice_color,'FaceAlpha',0.3,'EdgeAlpha',0.0)
% 
%     % plot the nominal value
% %     plot3(f_nom(1,i),f_nom(2,i),f_nom(3,i),'xk')
% end
% 
% 
% fric_cone_filename = 'Friction_Cone_Animation_v1.gif';
% 
% for i = 1:length(t_steps)
% 
%     % plot the friction cone
%     r = linspace(0,1,10);
%     theta = linspace(0,2*pi,50);
%     [RR,Theta] = meshgrid(r,theta);
%     X = RR.*cos(Theta);
%     Y = RR.*sin(Theta);
%     Z = 1./u_s.*RR;
%     h1 = surf(X,Y,Z,'EdgeColor','none','FaceColor','r','FaceAlpha','0.05');
%     
%     % plotting of the friction cone at the z-level
%     theta_friction = linspace(0,2*pi,100);
%     r_friction = f_nom(3,i)*u_s;
%     plot3(r_friction*cos(theta_friction),r_friction*sin(theta_friction),f_nom(3,i)*ones(1,length(theta_friction)),'-r')
% 
%     % plotting continuous force trajectory
%     % time
%     t_cont = linspace((i-1)/40,i/40); % 1/40 for a single iteration
%     % desired trajectory
%     for j = 1:length(t_cont)
%         [q_cont_des(:,j), qd_cont_des(:,j), qdd_cont_des(:,j)] = desired_trajectory(P, q_0, qd_0, qdd_0, t_cont(j), kvec);
%         % rnea
%         [u_temp f_temp n_temp] = rnea(q_cont_des(:,j), qd_cont_des(:,j), qd_cont_des(:,j), qdd_cont_des(:,j), true, params.nominal);
%     %     tau_int{i} = tau_temp{10,1};
%         f_cont(:,j) = f_temp(:,10);
%         n_cont(:,j) = n_temp(:,10);
%     end
%     plot3(f_cont(1,:),f_cont(2,:),f_cont(3,:),'-k', 'LineWidth', 3)
% 
%     % plot formatting
%     xlabel('x-axis Tangential Force (N)')
%     ylabel('y-axis Tangential Force (N)')
%     zlabel('z-axis Normal Force (N)')
%     grid on
%     axis('square')
%     make_animation( fig1,i,fric_cone_filename )
% end

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
    
    
    %         link_zonotopes{i,1} = zonotope([link_sizes(1,i)/2, link_sizes(1,i)/2; ...
    %                                         link_sizes(2,i)/2, link_sizes(2,i)/2; ...
    %                                         link_sizes(3,i)/2, link_sizes(3,i)/2]);
    
    
    
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

function plotConcentricCircles(fig, inner_radius, outer_radius, fill_color, edge_color)
% PLOTCONCENTRICCIRCLES plots two concentric circles with specified properties.
%   plotConcentricCircles(fig, inner_radius, outer_radius, fill_color, edge_color)
%
%   Inputs:
%       fig:          Figure handle where the circles will be plotted.
%       inner_radius: Inner radius of the donut (hole).
%       outer_radius: Outer radius of the donut (edge).
%       fill_color:   Color to fill the donut (e.g., 'r' for red).
%       edge_color:   Color of the donut's edge (e.g., 'b' for blue).

    % Switch to the specified figure or create a new one if not provided
    if nargin < 1
        fig = figure;
    else
        figure(fig); 
    end
    
    hold on;
    
    rectangle('Position',[-outer_radius,-outer_radius,2*outer_radius,2*outer_radius],...
        'Curvature',[1,1],...
        'FaceColor',fill_color,...
        'EdgeColor',edge_color); % Set the face color of the outer circle

    % Plot the inner circle
    rectangle('Position',[-inner_radius,-inner_radius,2*inner_radius,2*inner_radius],...
        'Curvature',[1,1],...
        'FaceColor','w',...     % Set the face color of the inner circle to white (hole)
        'EdgeColor',edge_color); % Set the edge color of the inner circle

%     hold off;

%     % Set axis limits to better view the donut
%     axis equal;
%     xlim([-outer_radius-1, outer_radius+1]);
%     ylim([-outer_radius-1, outer_radius+1]);
% 
%     % Set labels and title
%     xlabel('X');
%     ylabel('Y');
%     title('Concentric Circles');
% 
%     % Display the grid for better visualization (optional)
%     grid on;
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