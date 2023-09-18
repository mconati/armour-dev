
clear all;
close all;
clc;
dbstop if error

%% for agent
agent_urdf = 'Kinova_Grasp_w_Tray.urdf';
add_uncertainty_to = 'none'; % choose 'all', 'link', or 'none'
links_with_uncertainty = {}; % if add_uncertainty_to = 'link', specify links here.
uncertain_mass_range = [0.97, 1.03];

%% robot params:
robot = importrobot(agent_urdf);
robot.DataFormat = 'col';
robot.Gravity = [0 0 -9.81];
params = load_robot_params(robot, ...
                           'add_uncertainty_to', add_uncertainty_to, ...
                           'links_with_uncertainty', links_with_uncertainty,...
                           'uncertain_mass_range', uncertain_mass_range);

%% symbolic forward kinematics expression

q_plan = [sym('q1');
    sym('q2');
    sym('q3');
    sym('q4');
    sym('q5');
    sym('q6');
    sym('q7')];

fk_result = simplify(forward_kinematics_symbolic(q_plan,params.true.T0,params.true.joint_axes));

%% Extracting the Position Vector

ee_pos = fk_result(1:3,4);

q_test = double([pi/3;-pi/4;pi/12;pi/12;pi/4;-pi/3;pi/2]);

fk_result_eval = vpa(subs(fk_result,q_plan,q_test));

ee_pos_eval = fk_result_eval(1:3,4);

ee_pos_sym = fk_result(1:3,4);

ee_cost_sym = sum(ee_pos_sym);

%% Checking Against Real Forward Kinematics

fk_num_result = forward_kinematics(q_test,params.true.T0,params.true.joint_axes);
ee_pos_real = fk_result(1:3,4);

test_equivalent = ee_pos_eval - ee_pos_real;

%% Calculating partial derivative as test

diff(fk_result,q_plan(1));

%% Calculating the symbolic gradient
var1 = params.true.T0;
var2 = params.true.joint_axes;

grad_cost_ee = simplify(gradient(ee_cost_sym,q_plan));

eval_grad = vpa(subs(grad_cost_ee,q_plan,q_test))

grad_func = matlabFunction(grad_cost_ee);

ee_pos = fk_result(1:3,4);

cost_ee_sym = sum(ee_pos);

grad_cost_ee_sym = gradient(cost_sym)

