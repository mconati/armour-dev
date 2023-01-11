#include "Dynamics.h"

const std::string pathname = "/home/roahmlab/Documents/armour-dev/kinova_src/kinova_simulator_interfaces/kinova_planner_realtime_v2/buffer/";
const std::string inputfilename = pathname + "armour.in";
const std::string outputfilename1 = pathname + "armour.out";
const std::string outputfilename2 = pathname + "armour_joint_position_center.out";
const std::string outputfilename3 = pathname + "armour_joint_position_radius.out";
const std::string outputfilename4 = pathname + "armour_control_input_radius.out";
const std::string outputfilename5 = pathname + "armour_constraints.out";

int main() {
    /*
Section I:
    Parse input
    There is no check and warning, so be careful!
*/
    // Here is an example of required input
    // double q0[NUM_FACTORS] = {0.6543, -0.0876, -0.4837, -1.2278, -1.5735, -1.0720, 0};
    // double qd0[NUM_FACTORS] = {0, 0, 0, 0, 0, 0, 0};
    // double qdd0[NUM_FACTORS] = {0, 0, 0, 0, 0, 0, 0};
    // double q_des[NUM_FACTORS] = {0.6831, 0.009488, -0.2471, -0.9777, -1.414, -0.9958, 0};

    // const int num_obstacles = 10;
    // const double obstacles[num_obstacles * (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3] = {-0.28239,  -0.33281, 0.88069, 0.069825, 0, 0, 0,  0.09508, 0, 0, 0, 0.016624,
    //                                                                             -0.19033,  0.035391,  1.3032,  0.11024, 0, 0, 0, 0.025188, 0, 0, 0, 0.014342,
    //                                                                             0.67593, -0.085841, 0.43572,  0.17408, 0, 0, 0,  0.07951, 0, 0, 0,  0.18012,
    //                                                                             0.75382,   0.51895,  0.4731, 0.030969, 0, 0, 0,  0.22312, 0, 0, 0,  0.22981,
    //                                                                             0.75382,   0.51895,  0.4731, 0.030969, 0, 0, 0,  0.22312, 0, 0, 0,  0.22981,
    //                                                                             -0.28239,  -0.33281, 0.88069, 0.069825, 0, 0, 0,  0.09508, 0, 0, 0, 0.016624,
    //                                                                             -0.19033,  0.035391,  1.3032,  0.11024, 0, 0, 0, 0.025188, 0, 0, 0, 0.014342,
    //                                                                             0.67593, -0.085841, 0.43572,  0.17408, 0, 0, 0,  0.07951, 0, 0, 0,  0.18012,
    //                                                                             0.75382,   0.51895,  0.4731, 0.030969, 0, 0, 0,  0.22312, 0, 0, 0,  0.22981,
    //                                                                             0.75382,   0.51895,  0.4731, 0.030969, 0, 0, 0,  0.22312, 0, 0, 0,  0.22981};

    // declare this first and make sure we always have a new output
    std::ofstream outputstream1(outputfilename1);

    double q0[NUM_FACTORS] = {0.0};
    double qd0[NUM_FACTORS] = {0.0};
    double qdd0[NUM_FACTORS] = {0.0};
    double q_des[NUM_FACTORS] = {0.0};

    int num_obstacles = 0;
    double obstacles[MAX_OBSTACLE_NUM * (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3] = {0.0};

    std::ifstream inputstream(inputfilename);
    if (!inputstream.is_open()) {
        WARNING_PRINT("        CUDA & C++: Error reading input files !\n");
        outputstream1 << -1;
        outputstream1.close();
        throw;
    }
    for (int i = 0; i < NUM_FACTORS; i++) {
        inputstream >> q0[i];
    }
    for (int i = 0; i < NUM_FACTORS; i++) {
        inputstream >> qd0[i];
    }
    for (int i = 0; i < NUM_FACTORS; i++) {
        inputstream >> qdd0[i];
    }
    for (int i = 0; i < NUM_FACTORS; i++) {
        inputstream >> q_des[i];
    }
    inputstream >> num_obstacles;
    if (num_obstacles > MAX_OBSTACLE_NUM || num_obstacles < 0) {
        WARNING_PRINT("Number of obstacles larger than MAX_OBSTACLE_NUM !\n");
        outputstream1 << -1;
        outputstream1.close();
        throw;
    }
    if (num_obstacles > 0) {
        for (int i = 0; i < num_obstacles * (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3; i++) {
            inputstream >> obstacles[i];
        }
    }

    inputstream.close();

    double t_plan = 1.0; // optimize the distance between q_des and the desired trajectories at t_plan
     
/*
Section II:
    Initialize all polynomial zonotopes, including links and torques
*/
    double jointPositionRadius[NUM_TIME_STEPS * NUM_FACTORS * 3] = {0.0};

    // Obstacles O(obstacles, num_obstacles); 

    // do not count time for memory allocation
    auto start1 = std::chrono::high_resolution_clock::now();

    omp_set_num_threads(NUM_THREADS);
    int openmp_t_ind = 0; // openmp loop index

    /*
    Section II.A: Create JRS online
    */
    BezierCurve traj(q0, qd0, qdd0);

    try {
        // #pragma omp parallel for shared(traj) private(openmp_t_ind) schedule(static, NUM_TIME_STEPS / NUM_THREADS)
        for(openmp_t_ind = 0; openmp_t_ind < NUM_TIME_STEPS; openmp_t_ind++) {
            traj.makePolyZono(openmp_t_ind);
        }
    }
    catch (int errorCode) {
        WARNING_PRINT("Error creating JRS! Check previous error message!");
        return -1;
    }

    /*
    Section II.B: Compute link PZs and nominal torque PZs
    */
    KinematicsDynamics kd(&traj);

    try {
        // #pragma omp parallel for shared(kd) private(openmp_t_ind) schedule(static, NUM_TIME_STEPS / NUM_THREADS)
        for(openmp_t_ind = 0; openmp_t_ind < NUM_TIME_STEPS; openmp_t_ind++) {
            // compute nominal torque
            kd.rnea_nominal(openmp_t_ind);

            // compute interval torque
            kd.rnea_interval(openmp_t_ind);

            // compute max disturbance (stored in u_nom_int)
            for (int i = 0; i < NUM_FACTORS; i++) {
                kd.u_nom_int(i, openmp_t_ind) = kd.u_nom_int(i, openmp_t_ind) - kd.u_nom(i, openmp_t_ind);
            }

            // reduce non-only-k-dependent generators so that slice takes less time
            for (int i = 0; i < NUM_FACTORS; i++) {
                kd.u_nom(i, openmp_t_ind).reduce();
            }
        }
    }
    catch (int errorCode) {
        WARNING_PRINT("Error computing link PZs and nominal torque PZs! Check previous error message!");
        return -1;
    }

    /*
    Section II.C: Compute robust input bound
    */

    // the radius of the torque PZs
    Eigen::MatrixXd v_norm(NUM_FACTORS, NUM_TIME_STEPS);
    v_norm.setZero();

    try {
        for(int t_ind = 0; t_ind < NUM_TIME_STEPS; t_ind++) {
            // (1) add the bound of robust input
            Interval rho_max_temp = Interval(0.0);
            for (int i = 0; i < NUM_FACTORS; i++) {
                // compute norm of disturbance
                MatrixXInt temp = kd.u_nom_int(i, t_ind).toInterval(); // should be a 1-dim Interval
                rho_max_temp += temp(0) * temp(0);

                v_norm(i, t_ind) = alpha * (M_max - M_min) * eps + 0.5 * max(abs(temp(0).lower()), abs(temp(0).upper()));
            }
            rho_max_temp = sqrt(rho_max_temp);
            
            for (int i = 0; i < NUM_FACTORS; i++) {
                v_norm(i, t_ind) += 0.5 * rho_max_temp.upper();
            }

            // (2) add the radius of the nominal input PZ (after reducing)
            for (int i = 0; i < NUM_FACTORS; i++) {
                v_norm(i, t_ind) += kd.u_nom(i, t_ind).independent(0);
            }

            // (3) add friction
            for (int i = 0; i < NUM_FACTORS; i++) {
                v_norm(i, t_ind) += friction[i];
            }

            // so that v_norm would be the radius of the total control input PZ from now
        }
    }
    catch (int errorCode) {
        WARNING_PRINT("Error computing torque PZs! Check previous error message!");
        return -1;
    }

    // double factors[NUM_FACTORS] = {0.5, 0.6, 0.7, -0.5, -0.6, -0.7, 0.0};
    // for (int l = 0; l < NUM_FACTORS; l++) {
    //     MatrixXInt res = kd.u_nom(l, 63).slice(factors);
    //     cout << res;
    // }
}

// int main() {
//     Eigen::Matrix<double,1,3> center1;
//     center1 << 1,2,3;
//     PZsparse a(center1);

//     Eigen::Matrix<double,3,1> center2;
//     center2 << 4,5,6;
//     PZsparse b(center2);

//     cout << a << endl;

//     cout << b << endl;

//     PZsparse c = a * b;

//     // double factor[NUM_FACTORS] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7};
//     // MatrixXInt res = a.slice(factor);

//     // for (int i = 0; i < res.rows(); i++) {
//     //     for (int j = 0; j < res.cols(); j++) {
//     //         cout << "[ " << res(i,j).lower() << ", " << res(i,j).upper() << " ], ";
//     //     }
//     //     cout << endl;
//     // }

//     cout << c << endl;
    
//     return 0;
// }