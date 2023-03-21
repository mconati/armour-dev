%% Graph Forward Kinematics For Nodes
clear all;
close all;
clc;

%% Load Graph

graph_struct = load('QGraph700kPartial_AllInfo.mat');
Q_Graph = graph_struct.Q_Graph;

% Reduce Graph to Largest Connected Component

[bin, binsize] = conncomp(Q_Graph,'Type','weak');

idx = (binsize(bin) == max(binsize));

Q_Graph_Reduced = subgraph(Q_Graph,idx);

Q_Nodes = table2array(Q_Graph_Reduced.Nodes)';

%% Load Robot

agent_urdf = 'Kinova_Grasp_URDF.urdf';

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

%% Loop Through Joints and Calculate, Then Store Forward Kinematics

for i = 1:length(Q_Nodes)

    % calculate forward kinematics
    q_cur_fk = forward_kinematics_modified(Q_Nodes(:,i),params.true.T0,params.true.joint_axes);



end