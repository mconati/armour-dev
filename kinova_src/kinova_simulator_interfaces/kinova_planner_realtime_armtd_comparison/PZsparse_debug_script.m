clear; clc;

robot = importrobot('gen3.urdf');
robot.DataFormat = 'col';
robot.Gravity = [0 0 -9.81];
params = load_robot_params(robot, ...
                           'add_uncertainty_to', 'link', ...
                           'links_with_uncertainty', {'dumbbell_link'},...
                           'uncertain_mass_range', [0.97, 1.03]);

JRS_filenames = dir('../../src/offline_jrs/orig_parameterization/JRS_*.mat');
load('../../src/offline_jrs/orig_parameterization/c_kvi.mat');

q_0 = zeros(7,1);
q_dot_0 = linspace(-1,1,7);
q_des = zeros(7,1);

JRS_filename = {};
for i = 1:7
    [~, closest_idx] = min(abs(q_dot_0(i) - c_kvi));
    filename = sprintf('%sJRS_%0.3f.mat', '../../src/offline_jrs/orig_parameterization/', c_kvi(closest_idx));
    
    JRS_filename{end + 1} = filename;
end

cuda_input_file = fopen('results/armtd_main.in', 'w');

% q_0
for ind = 1:length(q_0)
    fprintf(cuda_input_file, '%.10f ', q_0(ind));
end
fprintf(cuda_input_file, '\n');

% q_dot_0
for ind = 1:length(q_dot_0)
    fprintf(cuda_input_file, '%.10f ', q_dot_0(ind));
end
fprintf(cuda_input_file, '\n');

% q_des
for ind = 1:length(q_des)
    fprintf(cuda_input_file, '%.10f ', q_des(ind));
end
fprintf(cuda_input_file, '\n');

k_range = zeros(7,1);

for ind = 1:7
    data = load(JRS_filename{ind});

    % cos(q)
    for i = 1:100
        c = data.JRS{i}.Z(1,1);
        fprintf(cuda_input_file, '%.10f ', c);
    end
    fprintf(cuda_input_file, '\n');
    for i = 1:100
        g = data.JRS{i}.Z(1,2);
        fprintf(cuda_input_file, '%.10f ', g);
    end
    fprintf(cuda_input_file, '\n');
    for i = 1:100
        r = sum(abs(data.JRS{i}.Z(1,3:end)), 2);
        fprintf(cuda_input_file, '%.10f ', r);
    end
    fprintf(cuda_input_file, '\n');
    
    % sin(q)
    for i = 1:100
        c = data.JRS{i}.Z(2,1);
        fprintf(cuda_input_file, '%.10f ', c);
    end
    fprintf(cuda_input_file, '\n');
    for i = 1:100
        g = data.JRS{i}.Z(2,2);
        fprintf(cuda_input_file, '%.10f ', g);
    end
    fprintf(cuda_input_file, '\n');
    for i = 1:100
        r = sum(abs(data.JRS{i}.Z(2,3:end)), 2);
        fprintf(cuda_input_file, '%.10f ', r);
    end
    fprintf(cuda_input_file, '\n');
    
    k_range(ind) = data.JRS{1}.Z(5,2);
    fprintf(cuda_input_file, '%.10f\n ', k_range(ind));
end

%%
for i = 1:7
    q_des(i) = orig_traj(1, q_0(i), q_dot_0(i), 0 * k_range(i));
end

JP = zeros(7,3);
for i = 1:7
    T = forward_kinematics(q_des(1:i), params.nominal.T0, params.nominal.joint_axes);
    JP(i,:) = T(1:3,4);
end

function [q_des, qd_des, qdd_des] = orig_traj(t, q_0, q_dot_0, k_scaled)

t_plan = 0.5;
t_stop = 1;

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

end
    
    