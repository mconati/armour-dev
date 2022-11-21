#include "NLPclass.h"
#include "BufferPath.h"

const std::string inputfilename = pathname + "armtd_main.in";
const std::string outputfilename1 = pathname + "armtd_main.out";
const std::string outputfilename2 = pathname + "armtd_main_joint_position_center.out";
const std::string outputfilename3 = pathname + "armtd_main_joint_position_radius.out";
const std::string outputfilename4 = pathname + "armtd_main_constraints.out";

int main() {
/*
Section I:
    Parse input
    There is no check and warning, so be careful!
*/
    // Here is an example of required input
    // TYPE q0[NUM_FACTORS] = {0.6543, -0.0876, -0.4837, -1.2278, -1.5735, -1.0720, 0};
    // TYPE qd0[NUM_FACTORS] = {0, 0, 0, 0, 0, 0, 0};
    // TYPE qdd0[NUM_FACTORS] = {0, 0, 0, 0, 0, 0, 0};
    // TYPE q_des[NUM_FACTORS] = {0.6831, 0.009488, -0.2471, -0.9777, -1.414, -0.9958, 0};

    // const int num_obstacles = 10;
    // const TYPE obstacles[num_obstacles * (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3] = {-0.28239,  -0.33281, 0.88069, 0.069825, 0, 0, 0,  0.09508, 0, 0, 0, 0.016624,
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

    TYPE q0[NUM_FACTORS] = {0.0};
    TYPE qd0[NUM_FACTORS] = {0.0};
    TYPE q_des[NUM_FACTORS] = {0.0};

    TYPE c_cos_q_des[NUM_FACTORS * NUM_TIME_STEPS];
    TYPE g_cos_q_des[NUM_FACTORS * NUM_TIME_STEPS];
    TYPE r_cos_q_des[NUM_FACTORS * NUM_TIME_STEPS];

    TYPE c_sin_q_des[NUM_FACTORS * NUM_TIME_STEPS];
    TYPE g_sin_q_des[NUM_FACTORS * NUM_TIME_STEPS];
    TYPE r_sin_q_des[NUM_FACTORS * NUM_TIME_STEPS];

    TYPE k_range[NUM_FACTORS];

    int num_obstacles = 0;
    TYPE obstacles[MAX_OBSTACLE_NUM * (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3] = {0.0};

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
        inputstream >> q_des[i];
    }
    for (int i = 0; i < NUM_FACTORS; i++) {
        for (int j = 0; j < NUM_TIME_STEPS; j++) {
            inputstream >> c_cos_q_des[i * NUM_TIME_STEPS + j];
        }
        for (int j = 0; j < NUM_TIME_STEPS; j++) {
            inputstream >> g_cos_q_des[i * NUM_TIME_STEPS + j];
        }
        for (int j = 0; j < NUM_TIME_STEPS; j++) {
            inputstream >> r_cos_q_des[i * NUM_TIME_STEPS + j];
        }
        for (int j = 0; j < NUM_TIME_STEPS; j++) {
            inputstream >> c_sin_q_des[i * NUM_TIME_STEPS + j];
        }
        for (int j = 0; j < NUM_TIME_STEPS; j++) {
            inputstream >> g_sin_q_des[i * NUM_TIME_STEPS + j];
        }
        for (int j = 0; j < NUM_TIME_STEPS; j++) {
            inputstream >> r_sin_q_des[i * NUM_TIME_STEPS + j];
        }
        inputstream >> k_range[i];
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
     
/*
Section II:
    Initialize all polynomial zonotopes, including forward kinematics and torques
*/
    PZsparse p[NUM_TIME_STEPS * NUM_FACTORS * 3];
    TYPE jointPositionRadius[NUM_TIME_STEPS * NUM_FACTORS * 3] = {0.0};

    ConstantAccelerationCurve traj(q0, qd0, 
                                   c_cos_q_des, g_cos_q_des, r_cos_q_des,
                                   c_sin_q_des, g_sin_q_des, r_sin_q_des,
                                   k_range);

    Obstacles O(obstacles, num_obstacles);  

    // do not count time for memory allocation
    auto start1 = std::chrono::high_resolution_clock::now();                                                         

    omp_set_num_threads(NUM_THREADS);

    int openmp_t_ind; // loop index

    #pragma omp parallel for shared(traj) private(openmp_t_ind) schedule(static, NUM_TIME_STEPS / NUM_THREADS)
    for(openmp_t_ind = 0; openmp_t_ind < NUM_TIME_STEPS; openmp_t_ind++) {
        // Part 1: convert parameterized Bezier curve desired trajectory to PZsparse
        traj.makePolyZono(openmp_t_ind);

        // Part 2: compute forward kinematics reachable sets
        KinematicsDynamics kd(traj.cos_q_des + openmp_t_ind * NUM_FACTORS, traj.sin_q_des + openmp_t_ind * NUM_FACTORS);

        kd.fk(p + openmp_t_ind * NUM_FACTORS * 3);

        for (int i = 0; i < NUM_FACTORS; i++) {
            p[(openmp_t_ind * NUM_FACTORS + i) * 3    ].reduce();
            p[(openmp_t_ind * NUM_FACTORS + i) * 3 + 1].reduce();
            p[(openmp_t_ind * NUM_FACTORS + i) * 3 + 2].reduce();

            jointPositionRadius[(openmp_t_ind * NUM_FACTORS + i) * 3    ] = getRadius(p[(openmp_t_ind * NUM_FACTORS + i) * 3    ].independent);
            jointPositionRadius[(openmp_t_ind * NUM_FACTORS + i) * 3 + 1] = getRadius(p[(openmp_t_ind * NUM_FACTORS + i) * 3 + 1].independent);
            jointPositionRadius[(openmp_t_ind * NUM_FACTORS + i) * 3 + 2] = getRadius(p[(openmp_t_ind * NUM_FACTORS + i) * 3 + 2].independent);
        }
    }

    O.initializeHyperPlane(jointPositionRadius);

    auto stop1 = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(stop1 - start1);
    cout << "        CUDA & C++: Time taken by generating trajectory & forward kinematics: " << duration1.count() << " milliseconds" << endl;

    // for (int j = 0; j < 3; j++) {
    //     p[((NUM_TIME_STEPS - 1) * NUM_FACTORS + 1) * 3 + j].print();
    // }

/*
Section III:
    Solve the optimization problem using IPOPT
*/
    auto start2 = std::chrono::high_resolution_clock::now();

    SmartPtr<armtd_NLP> mynlp = new armtd_NLP();
	mynlp->set_parameters(q_des, &traj, p, &O);

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-6);
	app->Options()->SetNumericValue("max_cpu_time", 10.0);
	app->Options()->SetIntegerValue("print_level", 0);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("linear_solver", "ma27");
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    // app->Options()->SetStringValue("output_file", "ipopt.out");
    // app->Options()->SetStringValue("derivative_test", "first-order");
    // app->Options()->SetNumericValue("derivative_test_perturbation", 1e-8);
    // app->Options()->SetNumericValue("derivative_test_tol", 1e-6);

    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if( status != Solve_Succeeded ) {
		WARNING_PRINT("Error during initialization!");
        outputstream1 << -1;
        outputstream1.close();
        throw;
    }

    // Ask Ipopt to solve the problem
	status = app->OptimizeTNLP(mynlp);

    if (status == Maximum_CpuTime_Exceeded) {
        cout << "        CUDA & C++: Ipopt maximum CPU time exceeded!\n";
    }

    auto stop2 = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(stop2 - start2);
    cout << "        CUDA & C++: Time taken by Ipopt: " << duration2.count() << " milliseconds" << endl;

/*
Section IV:
    Prepare output
*/
    // output k_opt
    // set precision to 10 decimal digits
    outputstream1 << std::setprecision(10);

    if (mynlp->feasible) {
        for (int i = 0; i < NUM_FACTORS; i++) {
            outputstream1 << mynlp->solution[i] << '\n';
        }
    }
    else {
        outputstream1 << -1 << '\n';
    }

    // output time cost (in milliseconds) in C++
    outputstream1 << duration1.count() + duration2.count();

    outputstream1.close();

    // output FRS and other information, you can comment them if they are unnecessary
    std::ofstream outputstream2(outputfilename2);
    outputstream2 << std::setprecision(10);
    for (int i = 0; i < NUM_TIME_STEPS * NUM_FACTORS; i++) {
        for (int j = 0; j < 3; j++) {
            outputstream2 << mynlp->checkJointsPosition[i * 3 + j] << ' ';
        }
        outputstream2 << '\n';
    }
    outputstream2.close();

    std::ofstream outputstream3(outputfilename3);
    outputstream3 << std::setprecision(6);
    for (int i = 0; i < NUM_TIME_STEPS * NUM_FACTORS; i++) {
        for (int j = 0; j < 3; j++) {
            outputstream3 << jointPositionRadius[i * 3 + j] << ' ';
        }
        outputstream3 << '\n';
    }
    outputstream3.close();

    std::ofstream outputstream4(outputfilename4);
    outputstream4 << std::setprecision(6);
    for (int i = 0; i < mynlp->constraint_number; i++) {
        outputstream4 << mynlp->g_copy[i] << '\n';
    }
    outputstream4.close();

    return 0;
}
