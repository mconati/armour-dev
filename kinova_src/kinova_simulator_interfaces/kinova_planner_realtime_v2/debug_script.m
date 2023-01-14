% clear; clc;
clc;

robot = importrobot('kinova_without_gripper.urdf');
robot.DataFormat = 'col';
robot.Gravity = [0 0 -9.81];
params = load_robot_params(robot, ...
                           'add_uncertainty_to', 'all', ...
                           'uncertain_mass_range', [0.97, 1.03]);

link_poly_zonotopes = create_pz_bounding_boxes(robot);

q0 = [0.3971798960 0.7189364630 0.4907726184 -0.9124350006 1.9214692866 0.7578017275 -0.3260810233]';
qd0 = [1.0000000000 1.0000000000 1.0000000000 1.0000000000 1.0000000000 1.0000000000 1.0000000000]';
qdd0 = [-1.0000000000 -1.0000000000 -1.0000000000 -1.0000000000 -1.0000000000 -1.0000000000 -1.0000000000]';

k = [0.5, 0.6, 0.7, -0.5, -0.6, -0.7, 0.0]';

k_range = [pi/24, pi/24, pi/24, pi/24, pi/24, pi/24, pi/30]';

q1 = q0 + k .* k_range;
qd1 = zeros(7,1);
qdd1 = zeros(7,1);

beta = match_deg5_bernstein_coefficients({q0, qd0, qdd0, q1, qd1, qdd1});

tspan = linspace(0, 1, 128 + 1);

for i = 128:128
    t_lb = tspan(i);
    t_ub = tspan(i + 1);

    t = (t_ub - t_lb) * rand + t_lb;

    [B, dB, ddB] = Bezier_kernel_deg5(t);

    q = zeros(7,1);
    qd = zeros(7,1);
    qdd = zeros(7,1);
    for j = 1:6
        q = q + beta{j} * B(j);
        qd = qd + beta{j} * dB(j);
        qdd = qdd + beta{j} * ddB(j);
    end

    [u, f, n] = rnea(q, qd, qd, qdd, true, params.nominal);
end