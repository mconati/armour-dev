

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

q_plan = [sym('q1'),...
    sym('q2'),...
    sym('q3'),...
    sym('q4'),...
    sym('q5'),...
    sym('q6'),...
    sym('q7')];

fk_result = simplify(forward_kinematics_symbolic(q_plan,params.true.T0,params.true.joint_axes));

%%

ee_pos = fk_result(1:3,4);

cost_ee_sym = sum(ee_pos);

grad_cost_ee_sym = gradient(cost_sym)