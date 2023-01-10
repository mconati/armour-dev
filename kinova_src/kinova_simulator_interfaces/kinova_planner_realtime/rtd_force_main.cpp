#include "NLPclass.h"
#include "BufferPath.h"

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
    TYPE qdd0[NUM_FACTORS] = {0.0};
    TYPE q_des[NUM_FACTORS] = {0.0};

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

    TYPE t_plan = 1.0; // optimize the distance between q_des and the desired trajectories at t_plan
     
/*
Section II:
    Initialize all polynomial zonotopes, including forward kinematics and torques
*/
    PZsparse p[NUM_TIME_STEPS * NUM_FACTORS * 3];
    PZsparse mass_arr[NUM_JOINTS];
    matPZsparse I_nominal_arr[NUM_JOINTS];
    matPZsparse I_uncertain_arr[NUM_JOINTS];
    PZsparse u_nom[NUM_TIME_STEPS * NUM_FACTORS];
    PZsparse u_nom_int[NUM_TIME_STEPS * NUM_FACTORS];
    PZsparse r[NUM_FACTORS];
    PZsparse Mr[NUM_TIME_STEPS * NUM_FACTORS];
    // PZsparse V[NUM_TIME_STEPS];
    TYPE jointPositionRadius[NUM_TIME_STEPS * NUM_FACTORS * 3] = {0.0};
    vecPZsparse f_c_nom[NUM_TIME_STEPS]; // may not be needed (can pass in nullptr instead)
    vecPZsparse n_c_nom[NUM_TIME_STEPS]; // may not be needed (can pass in nullptr instead)
    vecPZsparse f_c_int[NUM_TIME_STEPS];
    vecPZsparse n_c_int[NUM_TIME_STEPS];
    const Number u_s = 0.5; // static coefficient of friction between tray and object
    const Number surf_rad = 0.0762; // radius of contact area between tray and object (area assumed to be circular)
    // Note: might want to change this to be input to the C++ code from matlab?
    
    for (int i = 0; i < NUM_JOINTS; i++) {
        mass_arr[i] = PZsparse(mass[i], mass_uncertainty);
        I_nominal_arr[i] = matPZsparse(inertia + i * 9);
        I_uncertain_arr[i] = matPZsparse(inertia + i * 9, inertia_uncertainty);

        if (i < NUM_FACTORS) {
            r[i] = PZsparse(Interval(-eps, eps));
        }
    }

    Obstacles O(obstacles, num_obstacles); 

    // do not count time for memory allocation
    auto start1 = std::chrono::high_resolution_clock::now();

    BezierCurve traj(q0, qd0, qdd0);

    omp_set_num_threads(NUM_THREADS);

    int openmp_t_ind; // loop index

    // solving forward kinematics and dynamics in parallel
    #pragma omp parallel for shared(traj, mass, mass_arr, I_nominal_arr, I_uncertain_arr, r) private(openmp_t_ind) schedule(static, NUM_TIME_STEPS / NUM_THREADS)
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

        // Part 3: compute the nominal control input with nominal parameters
        kd.rnea(traj.qd_des + openmp_t_ind * NUM_FACTORS,
                traj.qda_des + openmp_t_ind * NUM_FACTORS,
                traj.qdda_des + openmp_t_ind * NUM_FACTORS,
                nullptr, // nominal mass already defined as constant array
                I_nominal_arr,
                u_nom + openmp_t_ind * NUM_FACTORS,
                nullptr,
                nullptr,
                true);
                // should output the nominal f_c and n_c and compare to matlab

        // Part 4: compute the nominal control input with interval parameters
        kd.rnea(traj.qd_des + openmp_t_ind * NUM_FACTORS,
                traj.qda_des + openmp_t_ind * NUM_FACTORS,
                traj.qdda_des + openmp_t_ind * NUM_FACTORS,
                mass_arr,
                I_uncertain_arr,
                u_nom_int + openmp_t_ind * NUM_FACTORS,
                f_c_int + openmp_t_ind, // ? Bohao question: not sure if the + openmp_t_ind is needed?
                n_c_int + openmp_t_ind, // ? Bohao question: not sure if the + openmp_t_ind is needed?  // don't need  * NUM_FACTORS since not storing for each joint
                true);

        // // Part 5: compute the Lyapunov function bound
        // // Part 5.a: compute M*r
        // kd.rnea(nullptr, // nullptr means zero vector
        //         nullptr, // nullptr means zero vector
        //         r,
        //         mass_arr,
        //         I_uncertain_arr,
        //         Mr + openmp_t_ind * NUM_FACTORS);
        //         // do not apply gravity here!
            
        // // Part 5.b: compute V = r^T * (M*r)
        // for (int i = 0; i < NUM_FACTORS; i++) {
        //     V[openmp_t_ind] += r[i] * Mr[openmp_t_ind * NUM_FACTORS + i];
        // }

        // // Part 5.c: compute V_diff = V - V
        // V[openmp_t_ind] = V[openmp_t_ind] - V[openmp_t_ind];

        // Part 6: compute rho_max
        for (int i = 0; i < NUM_FACTORS; i++) {
            // Part 6.a: compute disturbance
            u_nom_int[openmp_t_ind * NUM_FACTORS + i] = u_nom_int[openmp_t_ind * NUM_FACTORS + i] - u_nom[openmp_t_ind * NUM_FACTORS + i];
        }

        // Part 9: reduce non-only-k-dependent generators so that slice takes less time
        for (int i = 0; i < NUM_FACTORS; i++) {
            u_nom[openmp_t_ind * NUM_FACTORS + i].reduce();
        }

        // force constraints
        // Bohao question: need to reduce these?

    }

    O.initializeHyperPlane(jointPositionRadius);

    TYPE rho_max[NUM_TIME_STEPS] = {0.0};
    TYPE v_norm[NUM_TIME_STEPS * NUM_FACTORS] = {0.0};
   
    for(int t_ind = 0; t_ind < NUM_TIME_STEPS; t_ind++) {
        // // Part 5.c: compute V_diff = V - V
        // Interval V_diff_int = V[t_ind].toInterval();

        // // Part 6: compute rho_max
        // // Part 7: compute robust input bound
        // // ||v|| = alpha * ||V_diff|| / eps + 0.5 * R_i + ||R||
        Interval rho_max_temp = Interval(0.0);
        for (int i = 0; i < NUM_FACTORS; i++) {
            // Part 6.b: compute norm of disturbance
            Interval temp = u_nom_int[t_ind * NUM_FACTORS + i].toInterval();
            rho_max_temp += temp * temp;

            v_norm[t_ind * NUM_FACTORS + i] = alpha * (M_max - M_min) * eps + 0.5 * max(abs(temp.lower()), abs(temp.upper()));
        }
        rho_max_temp = sqrt(rho_max_temp);
        
        for (int i = 0; i < NUM_FACTORS; i++) {
            v_norm[t_ind * NUM_FACTORS + i] = v_norm[t_ind * NUM_FACTORS + i] + 0.5 * rho_max_temp.upper();
        }

        // Part 8: add the radius of the nominal input PZ, 
        //         so v_norm would be the radius of the total control input PZ from now
        for (int i = 0; i < NUM_FACTORS; i++) {
            v_norm[t_ind * NUM_FACTORS + i] = v_norm[t_ind * NUM_FACTORS + i] + getRadius(u_nom[t_ind * NUM_FACTORS + i].independent);
        }
    }

    auto stop1 = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(stop1 - start1);
    cout << "        CUDA & C++: Time taken by generating trajectory & rnea: " << duration1.count() << " milliseconds" << endl;

/*
Section III:
    Solve the optimization problem using IPOPT
*/
    auto start2 = std::chrono::high_resolution_clock::now();

    SmartPtr<armtd_NLP> mynlp = new armtd_NLP();
	mynlp->set_parameters(q_des, t_plan, &traj, p, u_nom, v_norm, &O, f_c_int, n_c_int, u_s, surf_rad);

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", IPOPT_OPTIMIZATION_TOLERANCE);
	app->Options()->SetNumericValue("max_cpu_time", IPOPT_MAX_CPU_TIME);
	app->Options()->SetIntegerValue("print_level", IPOPT_PRINT_LEVEL);
    app->Options()->SetStringValue("mu_strategy", IPOPT_MU_STRATEGY);
    app->Options()->SetStringValue("linear_solver", IPOPT_LINEAR_SOLVER);
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    // For gradient checking
    // app->Options()->SetStringValue("output_file", "ipopt.out");
    // app->Options()->SetStringValue("derivative_test", "first-order");
    // app->Options()->SetNumericValue("derivative_test_perturbation", 1e-8);
    // app->Options()->SetNumericValue("derivative_test_tol", 1e-6);

    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if( status != Solve_Succeeded ) {
		WARNING_PRINT("Error during initialization!");
        outputstream1 << -1 << '\n';
        outputstream1.close();
        throw;
    }

    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(mynlp);
	
    auto stop2 = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(stop2 - start2);

    if (status == Maximum_CpuTime_Exceeded) {
        cout << "        CUDA & C++: Ipopt maximum CPU time exceeded!\n";
    }
    
    if (status == Invalid_Option) {
        cout << "        CUDA & C++: Cannot find HSL library! Need to put libcoinhsl.so in proper path!\n";
    }
    else {
        cout << "        CUDA & C++: Time taken by Ipopt: " << duration2.count() << " milliseconds" << endl;
    }

/*
Section IV:
    Prepare output
*/
    // set precision to 10 decimal digits
    outputstream1 << std::setprecision(10);

    // output k_opt
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
    // std::ofstream outputstream2(outputfilename2);
    // outputstream2 << std::setprecision(10);
    // for (int i = 0; i < NUM_TIME_STEPS * NUM_FACTORS; i++) {
    //     for (int j = 0; j < 3; j++) {
    //         outputstream2 << mynlp->checkJointsPosition[i * 3 + j] << ' ';
    //     }
    //     outputstream2 << '\n';
    // }
    // outputstream2.close();

    // std::ofstream outputstream3(outputfilename3);
    // outputstream3 << std::setprecision(6);
    // for (int i = 0; i < NUM_TIME_STEPS * NUM_FACTORS; i++) {
    //     for (int j = 0; j < 3; j++) {
    //         outputstream3 << jointPositionRadius[i * 3 + j] << ' ';
    //     }
    //     outputstream3 << '\n';
    // }
    // outputstream3.close();

    // std::ofstream outputstream4(outputfilename4);
    // outputstream4 << std::setprecision(6);
    // for (int i = 0; i < NUM_TIME_STEPS * NUM_FACTORS; i++) {
    //     outputstream4 << v_norm[i] << '\n';
    // }
    // outputstream4.close();

    std::ofstream outputstream5(outputfilename5);
    outputstream5 << std::setprecision(6);
    for (int i = 0; i < mynlp->constraint_number; i++) {
        outputstream5 << mynlp->g_copy[i] << '\n';
    }
    outputstream5.close();

    return 0;
}
