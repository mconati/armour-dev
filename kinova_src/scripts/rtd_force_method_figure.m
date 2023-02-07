clear; close all; clc;

%% setup 2d3link robot
robot_name = 'Kinova_Grasp_Cylinder_Edge.urdf';
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
LLC_info.ultimate_bound = 0.191;
LLC_info.Kr = 10;

% create link poly zonotopes
[link_poly_zonotopes, link_sizes, meshes] = create_link_poly_zonos(robot);
zono_order = 10;

%% setup plotting flags
plot_idx = 0;
plot_trajectory_1 = true;
plot_trajectory_2 = false;
plot_trajectory_3 = false;
plot_trajectory_4 = false;
plot_forward_occupancy = true;
plot_pz_time = false;

plot_force_trajectory = true;

save_plot = false;

time_color = 1/256*[207,207,196];

traj_color = 1/256*[146,197,222];
traj_err_color = 1/256*[246,10,22];

pz_err_color = 1/256*[255, 165, 0];

slice_step_color = 1/256*[255, 0, 0];

unsliced_color = 1/256*[90,200,243];
unsliced_w_err_color = 1/256*[72,181,163];
sliced_color = 1/256*[165,137,193];

blues = 1/256*[222,235,247;
158,202,225;
49,130,189];
color = blues(1, :);
slicecolor = blues(2, :);
slicehardcolor = 1/256*[238, 244, 250];
linkcolor = blues(3, :);
% goalcolor = 1/256*[218,165,32];
goalcolor = 1/256*[0,187,51];

time_to_slice = 6;

%% compute trajectories
% initial conditions
q_0 = [0;-pi/2;0;0;0;0;0];
qd_0= [0;0;0;0;0;0;0];
qdd_0 = [0;0;0;0;0;0;0];

% trajectory parameter
kvec = [0.6; -0.8; 0.5; -0.2; -0.4; 0.35; 0.34];

% create pz trajectories
joint_axes = [zeros(2, length(q_0)); ones(1, length(q_0))]; % Question: what is this?
taylor_degree = 5;
traj_type = 'bernstein';
add_ultimate_bound = true;

[Q_des, Qd_des, Qdd_des, Q, Qd, Qd_a, Qdd_a, R_des, R_t_des, R, R_t, T, E_p, jrs_info] = create_jrs_online_modified(q_0,qd_0,qdd_0, joint_axes, taylor_degree, traj_type, add_ultimate_bound, LLC_info);

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
t_steps = 0:0.01:1;
n_steps = length(t_steps_dt);
q_des = zeros(length(kvec), length(t_steps));
q_des_dt = zeros(length(kvec), n_steps);
qd_des = zeros(length(kvec), n_steps);
qdd_des = zeros(length(kvec), n_steps);

for i = 1:101
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

for i = 1:jrs_info.n_t
    [tau_temp, f_temp, n_temp] = rnea(q_des, qd_des, qdd_des, true, params.nominal);
%     tau_int{i} = tau_temp{10,1};
    f_nom{i} = f_temp{10,1};
    n_nom{i} = n_temp{10,1};
end

%% Plotting Wrench Trajectory

if plot_force_trajectory

    plot_idx = plot_idx + 1;
    figure(plot_idx); clf; hold on;
    title('Force Plot: Last Joint')

    for i = 1:length(t_traj)
        % plot the polynomial overapproximation
        % calculate the inf/sup
        if f_int{1,i}.G
            poly_inf = f_int{1,i}.c(1) - sum(abs(f_int{1,i}.G(1,:))) - sum(abs(f_int{1,i}.Grest(1,:)));
            poly_sup = f_int{1,i}.c(1) + sum(abs(f_int{1,i}.G(1,:))) + sum(abs(f_int{1,i}.Grest(1,:)));
        else % the magnitude of Grest means that only those were tracked
            poly_inf = f_int{1,i}.c(1) - sum(abs(f_int{1,i}.Grest(1,:)));
            poly_sup = f_int{1,i}.c(1) + sum(abs(f_int{1,i}.Grest(1,:)));
        end
        p1 = patch([t_traj(i)+jrs_info.dt; t_traj(i)+jrs_info.dt; t_traj(i); t_traj(i)], [poly_sup; poly_inf; poly_inf; poly_sup], 'b');
%         p1.EdgeColor = pz_err_color;
        p1.LineWidth = 0.1;
%         p1.FaceColor = pz_err_color;

        % plot the nominal values
        
    end

end

%% Plotting Friction Cone

%% Plotting ZMP Diagram

%% Plotting Separation Constraint

% could maybe plot a 1D interval along the z-axis of the tray?

%% Notes
% Plotting the 3D force PZ on the robot tray would give info about both the
% separation and slipping constraints. Would be really cool visual.

%% plot trajectory 1
if plot_trajectory_1
    for j = 1:length(kvec)
        plot_idx = plot_idx + 1;
        figure(plot_idx); clf; hold on;

        % plot pz trajectories
        for i = 1:length(t_traj)
            % plot error polynomial zonotope interval
            poly_inf = Q_e{i, 1}{j, 1}.c - sum(abs(Q_e{i, 1}{j, 1}.G)) - sum(abs(Q_e{i, 1}{j, 1}.Grest));
            poly_sup = Q_e{i, 1}{j, 1}.c + sum(abs(Q_e{i, 1}{j, 1}.G)) + sum(abs(Q_e{i, 1}{j, 1}.Grest));
            p1 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_sup; poly_inf; poly_inf; poly_sup], 'b');
            p1.EdgeColor = pz_err_color;
            p1.LineWidth = 0.1;
            p1.FaceColor = pz_err_color;

            % plot unslice polynomial zonotope interval
            poly_inf = Q_des{i, 1}{j, 1}.c - sum(abs(Q_des{i, 1}{j, 1}.G)) - sum(abs(Q_des{i, 1}{j, 1}.Grest));
            poly_sup = Q_des{i, 1}{j, 1}.c + sum(abs(Q_des{i, 1}{j, 1}.G)) + sum(abs(Q_des{i, 1}{j, 1}.Grest));
            p2 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_sup; poly_inf; poly_inf; poly_sup], 'b');
            p2.EdgeColor = slice_step_color;
            p2.LineWidth = 0.1;
            p2.FaceColor = slice_step_color;

            % plot unsliced (original + err) polynomial zonotope interval
            poly_inf = Q{i, 1}{j, 1}.c - sum(abs(Q{i, 1}{j, 1}.G)) - sum(abs(Q{i, 1}{j, 1}.Grest));
            poly_sup = Q{i, 1}{j, 1}.c + sum(abs(Q{i, 1}{j, 1}.G)) + sum(abs(Q{i, 1}{j, 1}.Grest));
            p3 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_sup; poly_inf; poly_inf; poly_sup], 'b');
            p3.EdgeColor = slice_step_color;
            p3.LineWidth = 0.1;
            p3.FaceColor = slice_step_color;

            % plot sliced (original + err) polynomial zonotope interval
            poly_slice = getSubset(Q{i, 1}{j, 1}, id_slice(j), kvec(j));
            poly_slice_inf = poly_slice.c - sum(abs(poly_slice.G)) - sum(abs(poly_slice.Grest));
            poly_slice_sup = poly_slice.c + sum(abs(poly_slice.G)) + sum(abs(poly_slice.Grest));
            p4 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_slice_sup; poly_slice_inf; poly_slice_inf; poly_slice_sup], 'b');
            p4.EdgeColor = slice_step_color;
            p4.LineWidth = 0.1;
            p4.FaceColor = slice_step_color;

            % capture slice of time
            if i ~= time_to_slice
                p1.FaceAlpha = 0.2;
            end

        end

        axis tight;
        if plot_pz_time
            xl = xlim;
            yl = ylim-0.02;
            for i = 1:length(t_traj)
                
                poly_inf = yl(1)*T{end, 1}.c - sum(abs(T{end, 1}.G)) - sum(abs(T{end, 1}.Grest));
                poly_sup = yl(1)*T{end, 1}.c + sum(abs(T{end, 1}.G)) + sum(abs(T{end, 1}.Grest));
                p0 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_sup; poly_inf; poly_inf; poly_sup], 'b');
                p0.EdgeColor = time_color;
                p0.LineWidth = 0.1;
                p0.FaceColor = time_color;
            end
        
            axis tight;
            hold off;
        end
                % plot scalar trajectories
        plot(t_steps, q_des(j,:));
        plot(t_steps, q_max(j,:), 'Color', 'r');
        plot(t_steps, q_min(j,:), 'Color', 'r');

        if save_plot
            filename = sprintf('/Users/kronos/Research/armour/traj_error_%i.eps',j);
    
            exportgraphics(gcf,filename,'ContentType','vector','BackgroundColor','none');
        end
    end
end

%% plot trajectory 2
if plot_trajectory_2
    for j = 1:length(kvec)
        plot_idx = plot_idx + 1;
        figure(plot_idx); clf; hold on;

        % plot pz trajectories
        for i = 1:length(t_traj)
            % plot error polynomial zonotope interval
            poly_inf = Q_e{i, 1}{j, 1}.c - sum(abs(Q_e{i, 1}{j, 1}.G)) - sum(abs(Q_e{i, 1}{j, 1}.Grest));
            poly_sup = Q_e{i, 1}{j, 1}.c + sum(abs(Q_e{i, 1}{j, 1}.G)) + sum(abs(Q_e{i, 1}{j, 1}.Grest));
            p1 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_sup; poly_inf; poly_inf; poly_sup], 'b');
            p1.EdgeColor = slice_step_color;
            p1.LineWidth = 0.1;
            p1.FaceColor = slice_step_color;

            % plot unslice polynomial zonotope interval
            poly_inf = Q_des{i, 1}{j, 1}.c - sum(abs(Q_des{i, 1}{j, 1}.G)) - sum(abs(Q_des{i, 1}{j, 1}.Grest));
            poly_sup = Q_des{i, 1}{j, 1}.c + sum(abs(Q_des{i, 1}{j, 1}.G)) + sum(abs(Q_des{i, 1}{j, 1}.Grest));
            p2 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_sup; poly_inf; poly_inf; poly_sup], 'b');
            p2.EdgeColor = unsliced_color;
            p2.LineWidth = 0.1;
            p2.FaceColor = unsliced_color;

            % plot unsliced (original + err) polynomial zonotope interval
            poly_inf = Q{i, 1}{j, 1}.c - sum(abs(Q{i, 1}{j, 1}.G)) - sum(abs(Q{i, 1}{j, 1}.Grest));
            poly_sup = Q{i, 1}{j, 1}.c + sum(abs(Q{i, 1}{j, 1}.G)) + sum(abs(Q{i, 1}{j, 1}.Grest));
            p3 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_sup; poly_inf; poly_inf; poly_sup], 'b');
            p3.EdgeColor = slice_step_color;
            p3.LineWidth = 0.1;
            p3.FaceColor = slice_step_color;

            % plot sliced (original + err) polynomial zonotope interval
            poly_slice = getSubset(Q{i, 1}{j, 1}, id_slice(j), kvec(j));
            poly_slice_inf = poly_slice.c - sum(abs(poly_slice.G)) - sum(abs(poly_slice.Grest));
            poly_slice_sup = poly_slice.c + sum(abs(poly_slice.G)) + sum(abs(poly_slice.Grest));
            p4 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_slice_sup; poly_slice_inf; poly_slice_inf; poly_slice_sup], 'b');
            p4.EdgeColor = slice_step_color;
            p4.LineWidth = 0.1;
            p4.FaceColor = slice_step_color;

            % capture slice of time
            if i ~= time_to_slice
                p2.FaceAlpha = 0.2;
            end

        end

        axis tight;
        if plot_pz_time
            xl = xlim;
            yl = ylim-0.02;
            for i = 1:length(t_traj)
                
                poly_inf = yl(1)*T{end, 1}.c - sum(abs(T{end, 1}.G)) - sum(abs(T{end, 1}.Grest));
                poly_sup = yl(1)*T{end, 1}.c + sum(abs(T{end, 1}.G)) + sum(abs(T{end, 1}.Grest));
                p0 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_sup; poly_inf; poly_inf; poly_sup], 'b');
                p0.EdgeColor = time_color;
                p0.LineWidth = 0.1;
                p0.FaceColor = time_color;
            end
        
            axis tight;
            hold off;
        end
                % plot scalar trajectories
        plot(t_steps, q_des(j,:));
        plot(t_steps, q_max(j,:), 'Color', 'r');
        plot(t_steps, q_min(j,:), 'Color', 'r');

        if save_plot
            filename = sprintf('/Users/kronos/Research/armour/traj_unsliced_%i.eps',j);
    
            exportgraphics(gcf,filename,'ContentType','vector','BackgroundColor','none');
        end
    end
end

%% plot trajector 3
if plot_trajectory_3
    for j = 1:length(kvec)
        plot_idx = plot_idx + 1;
        figure(plot_idx); clf; hold on;

        % plot pz trajectories
        for i = 1:length(t_traj)
            % plot error polynomial zonotope interval
            poly_inf = Q_e{i, 1}{j, 1}.c - sum(abs(Q_e{i, 1}{j, 1}.G)) - sum(abs(Q_e{i, 1}{j, 1}.Grest));
            poly_sup = Q_e{i, 1}{j, 1}.c + sum(abs(Q_e{i, 1}{j, 1}.G)) + sum(abs(Q_e{i, 1}{j, 1}.Grest));
            p1 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_sup; poly_inf; poly_inf; poly_sup], 'b');
            p1.EdgeColor = slice_step_color;
            p1.LineWidth = 0.1;
            p1.FaceColor = slice_step_color;

            % plot unslice polynomial zonotope interval
            poly_inf = Q_des{i, 1}{j, 1}.c - sum(abs(Q_des{i, 1}{j, 1}.G)) - sum(abs(Q_des{i, 1}{j, 1}.Grest));
            poly_sup = Q_des{i, 1}{j, 1}.c + sum(abs(Q_des{i, 1}{j, 1}.G)) + sum(abs(Q_des{i, 1}{j, 1}.Grest));
            p2 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_sup; poly_inf; poly_inf; poly_sup], 'b');
            p2.EdgeColor = slice_step_color;
            p2.LineWidth = 0.1;
            p2.FaceColor = slice_step_color;

            % plot unsliced (original + err) polynomial zonotope interval
            poly_inf = Q{i, 1}{j, 1}.c - sum(abs(Q{i, 1}{j, 1}.G)) - sum(abs(Q{i, 1}{j, 1}.Grest));
            poly_sup = Q{i, 1}{j, 1}.c + sum(abs(Q{i, 1}{j, 1}.G)) + sum(abs(Q{i, 1}{j, 1}.Grest));
            p3 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_sup; poly_inf; poly_inf; poly_sup], 'b');
            p3.EdgeColor = unsliced_w_err_color;
            p3.LineWidth = 0.1;
            p3.FaceColor = unsliced_w_err_color;

            % plot sliced (original + err) polynomial zonotope interval
            poly_slice = getSubset(Q{i, 1}{j, 1}, id_slice(j), kvec(j));
            poly_slice_inf = poly_slice.c - sum(abs(poly_slice.G)) - sum(abs(poly_slice.Grest));
            poly_slice_sup = poly_slice.c + sum(abs(poly_slice.G)) + sum(abs(poly_slice.Grest));
            p4 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_slice_sup; poly_slice_inf; poly_slice_inf; poly_slice_sup], 'b');
            p4.EdgeColor = slice_step_color;
            p4.LineWidth = 0.1;
            p4.FaceColor = slice_step_color;

            % capture slice of time
            if i ~= time_to_slice
                p3.FaceAlpha = 0.2;
            end

        end

        axis tight;
        if plot_pz_time
            xl = xlim;
            yl = ylim-0.02;
            for i = 1:length(t_traj)
                
                poly_inf = yl(1)*T{end, 1}.c - sum(abs(T{end, 1}.G)) - sum(abs(T{end, 1}.Grest));
                poly_sup = yl(1)*T{end, 1}.c + sum(abs(T{end, 1}.G)) + sum(abs(T{end, 1}.Grest));
                p0 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_sup; poly_inf; poly_inf; poly_sup], 'b');
                p0.EdgeColor = time_color;
                p0.LineWidth = 0.1;
                p0.FaceColor = time_color;
            end
        
            axis tight;
            hold off;
        end
                % plot scalar trajectories
        plot(t_steps, q_des(j,:));
        plot(t_steps, q_max(j,:), 'Color', 'r');
        plot(t_steps, q_min(j,:), 'Color', 'r');

        if save_plot
            filename = sprintf('/Users/kronos/Research/armour/traj_total_error_%i.eps',j);
    
            exportgraphics(gcf,filename,'ContentType','vector','BackgroundColor','none');
        end
    end
end

%% plot trajector 4
if plot_trajectory_4
    for j = 1:length(kvec)
        plot_idx = plot_idx + 1;
        figure(plot_idx); clf; hold on;

        % plot pz trajectories
        for i = 1:length(t_traj)
            % plot error polynomial zonotope interval
            poly_inf = Q_e{i, 1}{j, 1}.c - sum(abs(Q_e{i, 1}{j, 1}.G)) - sum(abs(Q_e{i, 1}{j, 1}.Grest));
            poly_sup = Q_e{i, 1}{j, 1}.c + sum(abs(Q_e{i, 1}{j, 1}.G)) + sum(abs(Q_e{i, 1}{j, 1}.Grest));
            p1 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_sup; poly_inf; poly_inf; poly_sup], 'b');
            p1.EdgeColor = slice_step_color;
            p1.LineWidth = 0.1;
            p1.FaceColor = slice_step_color;

            % plot unslice polynomial zonotope interval
            poly_inf = Q_des{i, 1}{j, 1}.c - sum(abs(Q_des{i, 1}{j, 1}.G)) - sum(abs(Q_des{i, 1}{j, 1}.Grest));
            poly_sup = Q_des{i, 1}{j, 1}.c + sum(abs(Q_des{i, 1}{j, 1}.G)) + sum(abs(Q_des{i, 1}{j, 1}.Grest));
            p2 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_sup; poly_inf; poly_inf; poly_sup], 'b');
            p2.EdgeColor = slice_step_color;
            p2.LineWidth = 0.1;
            p2.FaceColor = slice_step_color;

            % plot unsliced (original + err) polynomial zonotope interval
            poly_inf = Q{i, 1}{j, 1}.c - sum(abs(Q{i, 1}{j, 1}.G)) - sum(abs(Q{i, 1}{j, 1}.Grest));
            poly_sup = Q{i, 1}{j, 1}.c + sum(abs(Q{i, 1}{j, 1}.G)) + sum(abs(Q{i, 1}{j, 1}.Grest));
            p3 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_sup; poly_inf; poly_inf; poly_sup], 'b');
            p2.EdgeColor = slice_step_color;
            p3.LineWidth = 0.1;
            p3.FaceColor = slice_step_color;

            % plot sliced (original + err) polynomial zonotope interval
            poly_slice = getSubset(Q{i, 1}{j, 1}, id_slice(j), kvec(j));
            poly_slice_inf = poly_slice.c - sum(abs(poly_slice.G)) - sum(abs(poly_slice.Grest));
            poly_slice_sup = poly_slice.c + sum(abs(poly_slice.G)) + sum(abs(poly_slice.Grest));
            p4 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_slice_sup; poly_slice_inf; poly_slice_inf; poly_slice_sup], 'b');
            p4.EdgeColor = sliced_color;
            p4.LineWidth = 0.1;
            p4.FaceColor = sliced_color;

            % capture slice of time
            if i ~= time_to_slice
                p4.FaceAlpha = 0.2;
            end

        end

        axis tight;
        if plot_pz_time
            xl = xlim;
            yl = ylim-0.02;
            for i = 1:length(t_traj)
                
                poly_inf = yl(1)*T{end, 1}.c - sum(abs(T{end, 1}.G)) - sum(abs(T{end, 1}.Grest));
                poly_sup = yl(1)*T{end, 1}.c + sum(abs(T{end, 1}.G)) + sum(abs(T{end, 1}.Grest));
                p0 = patch([t_traj(i)+jrs_info.dt/2; t_traj(i)+jrs_info.dt/2; t_traj(i) - jrs_info.dt/2; t_traj(i) - jrs_info.dt/2], [poly_sup; poly_inf; poly_inf; poly_sup], 'b');
                p0.EdgeColor = time_color;
                p0.LineWidth = 0.1;
                p0.FaceColor = time_color;
            end
        
            axis tight;
            hold off;
        end
                % plot scalar trajectories
        plot(t_steps, q_des(j,:));
        plot(t_steps, q_max(j,:), 'Color', 'r');
        plot(t_steps, q_min(j,:), 'Color', 'r');

        if save_plot
            filename = sprintf('/Users/kronos/Research/armour/traj_slice_%i.eps',j);
    
            exportgraphics(gcf,filename,'ContentType','vector','BackgroundColor','none');
        end
    end
end

%% compute forward occupancy
for i = 1:jrs_info.n_t
    % desired
    [R_w_des{i, 1}, p_w_des{i, 1}] = pzfk(R_des{i, 1}, params.pz_nominal); 
    
    % actual
    [R_w{i, 1}, p_w{i, 1}] = pzfk(R{i, 1}, params.pz_nominal); 

    % error
    [R_e_w{i, 1}, p_e_w{i, 1}] = pzfk(R_e{i, 1}, params.pz_nominal); 

   for j = 1:length(link_poly_zonotopes)
      % unsliced desired 
      FO_des{i, 1}{j, 1} = R_w_des{i, 1}{j, 1}*link_poly_zonotopes{j, 1} + p_w_des{i, 1}{j, 1};
      FO_des{i, 1}{j, 1} = reduce(FO_des{i, 1}{j, 1}, 'girard', params.pz_interval.zono_order);
      FO_des{i, 1}{j, 1} = remove_dependence(FO_des{i, 1}{j, 1}, jrs_info.k_id(end));

      % sliced desired
      fully_sliceable_des_tmp = polyZonotope_ROAHM(FO_des{i}{j}.c, FO_des{i}{j}.G, [], FO_des{i}{j}.expMat, FO_des{i}{j}.id);
      FO_des_sliced{i,1}{j,1} = zonotope([slice(fully_sliceable_des_tmp, kvec), FO_des{i}{j}.Grest]);

      % unsliced actual
      FO{i, 1}{j, 1} = R_w{i, 1}{j, 1}*link_poly_zonotopes{j, 1} + p_w{i, 1}{j, 1}; 
      FO{i, 1}{j, 1} = reduce(FO{i, 1}{j, 1}, 'girard', params.pz_interval.zono_order);
      FO{i, 1}{j, 1} = remove_dependence(FO{i, 1}{j, 1}, jrs_info.k_id(end));

      % error
      FO_e{i, 1}{j, 1} = R_e_w{i, 1}{j, 1}*link_poly_zonotopes{j, 1} + p_e_w{i, 1}{j, 1}; 
      FO_e{i, 1}{j, 1} = reduce(FO_e{i, 1}{j, 1}, 'girard', params.pz_interval.zono_order);
      FO_e{i, 1}{j, 1} = remove_dependence(FO_e{i, 1}{j, 1}, jrs_info.k_id(end));

      % sliced actual
      fully_sliceable_tmp = polyZonotope_ROAHM(FO{i}{j}.c, FO{i}{j}.G, [], FO{i}{j}.expMat, FO{i}{j}.id);
      FO_sliced{i,1}{j,1} = zonotope([slice(fully_sliceable_tmp, kvec), FO{i}{j}.Grest]);
   end
end

%% plot forward occupancy 1
% plot time backward
if plot_forward_occupancy
    plot_idx = plot_idx + 1;
    figure(plot_idx); clf; hold on;
    for i = jrs_info.n_t:-1:1
         % plot forward occupancy
        for j = 1:length(link_poly_zonotopes)
            % pz ultimate bounds
            Z = zonotope(FO_e{i, 1}{j, 1});
            p1 = plot(Z, [1 2], 'Filled', true);
            p1.FaceColor = pz_err_color;

            % unsliced 
            Z = zonotope(FO_des{i, 1}{j, 1});
            p2 = plot(Z, [1 2], 'Filled', true);
            p2.FaceColor = slice_step_color;

            % unsliced with error
            Z = zonotope(FO{i, 1}{j, 1});
            p3 = plot(Z, [1 2], 'Filled', true);
            p3.FaceColor = slice_step_color;
            
            % sliced
            Z = zonotope(FO_sliced{i, 1}{j, 1});
            p4 = plot(Z, [1 2], 'Filled', true);
            p4.FaceColor = slice_step_color;

            % capture slice of time
            if i ~= time_to_slice
                p1.FaceAlpha = 0.2;
            end
        end

    end

    % setup links
    [xlinks, ylinks] = link_points(robot, q_des_dt(:,time_to_slice));
    [xee, yee] = ee_points(robot, q_des_dt(:,time_to_slice), link_sizes(1:2, end));

    % plot links
    plot(xlinks, ylinks, 'Color', linkcolor, 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 1);

    % plot end effector
    plot([xlinks(end),xee(1:2)], [ylinks(end),yee(1:2)], 'Color', linkcolor, 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 1);
    plot([xlinks(end),xee(3:4)], [ylinks(end),yee(3:4)], 'Color', linkcolor, 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 1);

    if save_plot
        filename = '/Users/kronos/Research/armour/fo_error.eps';
        exportgraphics(gcf,filename,'ContentType','vector','BackgroundColor','none');
    end
end

%% plot forward occupancy 2
% plot time backward
if plot_forward_occupancy
    plot_idx = plot_idx + 1;
    figure(plot_idx); clf; hold on;
    for i = jrs_info.n_t:-1:1
         % plot forward occupancy
        for j = 1:length(link_poly_zonotopes)
       
            % pz ultimate bounds
            Z = zonotope(FO_e{i, 1}{j, 1});
            p1 = plot(Z, [1 2], 'Filled', true);
            p1.FaceColor = slice_step_color;

            % unsliced 
            Z = zonotope(FO_des{i, 1}{j, 1});
            p2 = plot(Z, [1 2], 'Filled', true);
            p2.FaceColor = unsliced_color;

            % unsliced with error
            Z = zonotope(FO{i, 1}{j, 1});
            p3 = plot(Z, [1 2], 'Filled', true);
            p3.FaceColor = slice_step_color;
            
            % sliced
            Z = zonotope(FO_sliced{i, 1}{j, 1});
            p4 = plot(Z, [1 2], 'Filled', true);
            p4.FaceColor = slice_step_color;

            % capture slice of time
            if i ~= time_to_slice
                p2.FaceAlpha = 0.2;
            end
        end
    end

    % setup links
    [xlinks, ylinks] = link_points(robot, q_des_dt(:,time_to_slice));
    [xee, yee] = ee_points(robot, q_des_dt(:,time_to_slice), link_sizes(1:2, end));

    % plot links
    plot(xlinks, ylinks, 'Color', linkcolor, 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 1);

    % plot end effector
    plot([xlinks(end),xee(1:2)], [ylinks(end),yee(1:2)], 'Color', linkcolor, 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 1);
    plot([xlinks(end),xee(3:4)], [ylinks(end),yee(3:4)], 'Color', linkcolor, 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 1);

    if save_plot
        filename = '/Users/kronos/Research/armour/fo_nominal.eps';
        exportgraphics(gcf,filename,'ContentType','vector','BackgroundColor','none');
    end
end


%% plot forward occupancy 3
% plot time backward
if plot_forward_occupancy
    plot_idx = plot_idx + 1;
    figure(plot_idx); clf; hold on;
    for i = jrs_info.n_t:-1:1
         % plot forward occupancy
        for j = 1:length(link_poly_zonotopes)
       
            % pz ultimate bounds
            Z = zonotope(FO_e{i, 1}{j, 1});
            p1 = plot(Z, [1 2], 'Filled', true);
            p1.FaceColor = slice_step_color;

            % unsliced 
            Z = zonotope(FO_des{i, 1}{j, 1});
            p2 = plot(Z, [1 2], 'Filled', true);
            p2.FaceColor = slice_step_color;

            % unsliced with error
            Z = zonotope(FO{i, 1}{j, 1});
            p3 = plot(Z, [1 2], 'Filled', true);
            p3.FaceColor = unsliced_w_err_color;
            
            % sliced
            Z = zonotope(FO_sliced{i, 1}{j, 1});
            p4 = plot(Z, [1 2], 'Filled', true);
            p4.FaceColor = slice_step_color;

            % capture slice of time
            if i ~= time_to_slice
                p3.FaceAlpha = 0.2;
            end
        end
    end

    % setup links
    [xlinks, ylinks] = link_points(robot, q_des_dt(:,time_to_slice));
    [xee, yee] = ee_points(robot, q_des_dt(:,time_to_slice), link_sizes(1:2, end));

    % plot links
    plot(xlinks, ylinks, 'Color', linkcolor, 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 1);

    % plot end effector
    plot([xlinks(end),xee(1:2)], [ylinks(end),yee(1:2)], 'Color', linkcolor, 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 1);
    plot([xlinks(end),xee(3:4)], [ylinks(end),yee(3:4)], 'Color', linkcolor, 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 1);

    if save_plot
        filename = '/Users/kronos/Research/armour/fo_total_error.eps';
        exportgraphics(gcf,filename,'ContentType','vector','BackgroundColor','none');
    end
end


%% plot forward occupancy 4
% plot time backward
if plot_forward_occupancy
    plot_idx = plot_idx + 1;
    figure(plot_idx); clf; hold on;
    for i = jrs_info.n_t:-1:1
         % plot forward occupancy
        for j = 1:length(link_poly_zonotopes)
       
            % pz ultimate bounds
            Z = zonotope(FO_e{i, 1}{j, 1});
            p1 = plot(Z, [1 2], 'Filled', true);
            p1.FaceColor = slice_step_color;

            % unsliced 
            Z = zonotope(FO_des{i, 1}{j, 1});
            p2 = plot(Z, [1 2], 'Filled', true);
            p2.FaceColor = slice_step_color;

            % unsliced with error
            Z = zonotope(FO{i, 1}{j, 1});
            p3 = plot(Z, [1 2], 'Filled', true);
            p3.FaceColor = slice_step_color;
            
            % sliced
            Z = zonotope(FO_sliced{i, 1}{j, 1});
            p4 = plot(Z, [1 2], 'Filled', true);
            p4.FaceColor = sliced_color;

            % capture slice of time
            if i ~= time_to_slice
                p4.FaceAlpha = 0.2;
            end
        end
    end

    % setup links
    [xlinks, ylinks] = link_points(robot, q_des_dt(:,time_to_slice));
    [xee, yee] = ee_points(robot, q_des_dt(:,time_to_slice), link_sizes(1:2, end));

    % plot links
    plot(xlinks, ylinks, 'Color', linkcolor, 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 1);

    % plot end effector
    plot([xlinks(end),xee(1:2)], [ylinks(end),yee(1:2)], 'Color', linkcolor, 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 1);
    plot([xlinks(end),xee(3:4)], [ylinks(end),yee(3:4)], 'Color', linkcolor, 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 1);

    if save_plot
        filename = '/Users/kronos/Research/armour/fo_slice.eps';
        exportgraphics(gcf,filename,'ContentType','vector','BackgroundColor','none');
    end
end
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
