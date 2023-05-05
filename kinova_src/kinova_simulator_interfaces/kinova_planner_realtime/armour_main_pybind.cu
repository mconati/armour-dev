#include "NLPclass.h"
#include "FastCollisionChecking.h"
#include "ReachsetsPath.h"
#include <cmath>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <memory>

namespace py = pybind11;

#define NUM_EDGES 112826403 // 50221615
#define COLLISION_THRESHOLD -0.075 // -0.10
#define EDGE_THRESHOLD 0.45 // 1000

class pzsparse {
    private:

        std::string data_dir;
        const std::string ofile1 = "/armour.out";
        const std::string ofile2 = "/armour_joint_position_center.out";
        const std::string ofile3 = "/armour_joint_position_radius.out";

        // Eigen::VectorXd a ({1.5, 2.5, 3.5});
        // Eigen::MatrixXd b ({1.5, 2.5, 3.5});

        Eigen::VectorXd q0 = Eigen::VectorXd::Zero(NUM_FACTORS);
        Eigen::VectorXd qd0 = Eigen::VectorXd::Zero(NUM_FACTORS);
        Eigen::VectorXd qdd0 = Eigen::VectorXd::Zero(NUM_FACTORS);
        Eigen::VectorXd q_des = Eigen::VectorXd::Zero(NUM_FACTORS);

        // Eigen::VectorXd q0; //q0.setZero(NUM_FACTORS);
        // Eigen::MatrixXd qd0; //qd0.setZero();
        // Eigen::MatrixXd qdd0; //qdd0.setZero();
        // Eigen::MatrixXd q_des; //q_des.setZero();
 
        double t_plan = 0.5*DURATION; // 0.5; // 0.5; // 
        // Kinova Hardware Demo Values: u_s = 0.609382421; surf_rad =  0.058/2;
        double u_s = 0.6;  // 0.396674703; // 0.3358; // 0.609382421; // 0.5; // static coefficient of friction between tray and object
        double surf_rad =  0.058 / 2; // 0.0762; // RADIUS of contact area between tray and object (area assumed to be circular) 
        // Note: might want to change this to be input to the C++ code from matlab?

        int num_obstacles = 0;
        double obstacles[MAX_OBSTACLE_NUM * (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3] = {0.0};
        std::shared_ptr<Obstacles> O_ptr{nullptr};
        std::shared_ptr<SimplifiedObstacles> O_ptr2{nullptr};
        // Obstacles O(obstacles, num_obstacles); 
        // Obstacles O;

        BezierCurve traj;
        PZsparse p[NUM_TIME_STEPS * NUM_FACTORS * 3];
        PZsparse u_nom[NUM_TIME_STEPS * NUM_FACTORS];
        double v_norm[NUM_TIME_STEPS * NUM_FACTORS];
        double jointPositionRadius[NUM_TIME_STEPS * NUM_FACTORS * 3];

        void set_obstacles(py::array_t<double> obstacle_vec){
            auto obstacle_ = obstacle_vec.unchecked<2>();
            num_obstacles = obstacle_.shape(0);
            int row = 0;
            int col = 0;
            int obs_dim = (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3;
            cout << "Obstacles set to:" << endl;
            for (int i = 0; i < num_obstacles * obs_dim; i++) {
                row = i / obs_dim;
                col = i % obs_dim;
                obstacles[i] = obstacle_(row, col);
            }
            for (int i = 0; i <= row; i++){
                cout << "obstacle" << i+1 << ": [ ";
                for (int j = 0; j <= col; j++){
                    cout << obstacles[i*obs_dim + j] << ' ';
                }
                cout << "] \n";
            }
            cout << '\n' << endl;
            std::cout << "allocating obstacles..." << std::endl;
            O_ptr = std::make_unique<Obstacles>();
            O_ptr->initialize(obstacles, num_obstacles);
            // O.initialize(obstacles, num_obstacles);

            O_ptr2 = std::make_unique<SimplifiedObstacles>();
            O_ptr2->initialize(obstacles, num_obstacles);

            std::cout << "Obstacles allocated!" << std::endl;
        }

        void set_goal(Eigen::Ref<Eigen::VectorXd> qdes_vec){
            auto qdes_ = qdes_vec; //.unchecked<1>();
            // cout << "Goal set to: \n [";
            for (int i = 0; i < NUM_FACTORS; i++) {
                q_des[i] = qdes_(i);
                // cout << q_des[i] << ' ';
            }
            // cout << "]" << endl;
        }

        void set_state(Eigen::Ref<Eigen::VectorXd> q0_vec, Eigen::Ref<Eigen::VectorXd> qd0_vec, Eigen::Ref<Eigen::VectorXd> qdd0_vec){
            auto q0_ = q0_vec; //.unchecked<1>();
            auto qd0_ = qd0_vec; //.unchecked<1>();
            auto qdd0_ = qdd0_vec; //.unchecked<1>();
            // cout << "States set to:" << endl;
            // cout << "[q0,     qd0,     qdd0]" << endl;
            for (int i = 0; i < NUM_FACTORS; i++){
                q0[i] = q0_(i);
                qd0[i] = qd0_(i);
                qdd0[i] = qdd0_(i);
                // cout << q0[i] << ' ' << qd0[i] << ' ' << qdd0[i] << endl;
            }
        }

        // void write_reachset(SmartPtr<armtd_NLP> mynlp){
        //     cout << "Saving reachsets..." << endl;
        //     for (int link = 0; link < 7; link++){
        //         std::ofstream file(reachset+std::to_string(link+1)+".txt");
        //         file << std::setprecision(8);
        //         for (int point = 0; point < 8; point++){
        //             int axis1 = static_cast<int>(point / 4);
        //             int axis2 = static_cast<int>(point / 2);
        //             int axis3 = point;
        //             int axis[] = {axis1, axis2, axis3};
        //             for (int i = 0; i < NUM_TIME_STEPS; i++) {
        //                 for (int j = 0; j < 3; j++) {
        //                     file << mynlp->checkJointsPosition[(i * NUM_FACTORS + link) * 3 + j] + (link_radius[link][0]+jointPositionRadius[(i * NUM_FACTORS + link) * 3 + j])*std::pow(-1, axis[j])<< ',';
        //                 }
        //                 file << '\n';
        //             }
        //         }
        //         file.close();
        //     }
        //     // write the origin
        //     std::ofstream file(reachset+"0.txt");
        //     for (int point = 0; point<8; point++){
        //         int axis1 = static_cast<int>(point / 4);
        //         int axis2 = static_cast<int>(point / 2);
        //         int axis3 = point;
        //         int axis[] = {axis1, axis2, axis3};
        //         for (int j = 0; j < 3; j++) {
        //             file << (0.04)*std::pow(-1, axis[j])<< ',';
        //         }
        //         file << '\n';
        //     }
        // }

    public:
    
        pzsparse(py::array_t<double> obs_vec, const std::string &dir)
        : data_dir(dir), num_obstacles(0){
            
            set_obstacles(obs_vec);
        }

        ~pzsparse()
        {
        }

        const int getNumObstacles(){
            return num_obstacles;
        }
        
        py::array_t<double> optimize(Eigen::Ref<Eigen::VectorXd> q0_vec, Eigen::Ref<Eigen::VectorXd> qd0_vec, Eigen::Ref<Eigen::VectorXd> qdd0_vec, 
                                    Eigen::Ref<Eigen::VectorXd> qgoal_vec) {
            

            set_goal(qgoal_vec);
            set_state(q0_vec, qd0_vec, qdd0_vec);      

            auto start1 = std::chrono::high_resolution_clock::now();      

            // Create JRS online
            traj = BezierCurve(q0, qd0, qdd0);
            omp_set_num_threads(NUM_THREADS);
            int openmp_t_ind = 0; // openmp loop index

            try {
                #pragma omp parallel for shared(traj) private(openmp_t_ind) schedule(static, NUM_TIME_STEPS / NUM_THREADS)
                for(openmp_t_ind = 0; openmp_t_ind < NUM_TIME_STEPS; openmp_t_ind++) {
                    traj.makePolyZono(openmp_t_ind);
                }
            }
            catch (int errorCode) {
                WARNING_PRINT("        CUDA & C++: Error creating JRS! Check previous error message!");
                throw;
            }

            // Compute link PZs and nominal torque PZs
            KinematicsDynamics kd(&traj);
            Eigen::Matrix<double, 3, 3 + 3> link_independent_generators[NUM_TIME_STEPS * NUM_JOINTS];

            try {
                #pragma omp parallel for shared(kd, link_independent_generators) private(openmp_t_ind) schedule(static, NUM_TIME_STEPS / NUM_THREADS)
                for(openmp_t_ind = 0; openmp_t_ind < NUM_TIME_STEPS; openmp_t_ind++) {
                    // compute link PZs through forward kinematics
                    kd.fk(openmp_t_ind);

                    // reduce non-only-k-dependent generators so that slice takes less time
                    for (int i = 0; i < NUM_JOINTS; i++) {
                        link_independent_generators[openmp_t_ind * NUM_JOINTS + i] = kd.links(i, openmp_t_ind).reduce_link_PZ();
                    }

                    // compute nominal torque
                    kd.rnea_nominal(openmp_t_ind);

                    // compute interval torque
                    kd.rnea_interval(openmp_t_ind);

                    // compute max disturbance (stored in u_nom_int)
                    for (int i = 0; i < NUM_FACTORS; i++) {
                        kd.u_nom_int(i, openmp_t_ind) = kd.u_nom_int(i, openmp_t_ind) - kd.u_nom(i, openmp_t_ind);
                        // kd.u_nom_int(i, openmp_t_ind) = 0;
                    }

                    // reduce non-only-k-dependent generators so that slice takes less time
                    for (int i = 0; i < NUM_FACTORS; i++) {
                        kd.u_nom(i, openmp_t_ind).reduce();
                    }
                }
            }
            catch (int errorCode) {
                WARNING_PRINT("        CUDA & C++: Error computing link PZs and nominal torque PZs! Check previous error message!");
                throw;
            }

            // Compute robust input bound
            Eigen::MatrixXd torque_radius(NUM_FACTORS, NUM_TIME_STEPS);
            torque_radius.setZero();

            try {
                for(int t_ind = 0; t_ind < NUM_TIME_STEPS; t_ind++) {
                    // (1) add the bound of robust input (||v||)
                    Interval rho_max_temp = Interval(0.0);
                    for (int i = 0; i < NUM_FACTORS; i++) {
                        // compute norm of disturbance
                        MatrixXInt temp = kd.u_nom_int(i, t_ind).toInterval(); // should be a 1-dim Interval
                        rho_max_temp += temp(0) * temp(0);

                        torque_radius(i, t_ind) = alpha * (M_max - M_min) * eps + 0.5 * max(abs(temp(0).lower()), abs(temp(0).upper()));
                        // torque_radius(i, t_ind) = 0;
                    }
                    rho_max_temp = sqrt(rho_max_temp);
                    
                    for (int i = 0; i < NUM_FACTORS; i++) {
                        torque_radius(i, t_ind) += 0.5 * rho_max_temp.upper();
                        // torque_radius(i, t_ind) += 0;
                    }

                    // (2) add the radius of the nominal input PZ (after reducing)
                    for (int i = 0; i < NUM_FACTORS; i++) {
                        torque_radius(i, t_ind) += kd.u_nom(i, t_ind).independent(0);
                        // torque_radius(i, t_ind) += 0;
                    }

                    // (3) add friction
                    for (int i = 0; i < NUM_FACTORS; i++) {
                        torque_radius(i, t_ind) += friction[i];
                        // torque_radius(i, t_ind) += 0;
                    }

                    // so that torque_radius would be the radius of the total control input PZ from now
                }
            }
            catch (int errorCode) {
                WARNING_PRINT("        CUDA & C++: Error computing torque PZs! Check previous error message!");
                throw;
            }

            // Buffer obstacles and initialize collision checking hyperplanes
            try {
                O_ptr->initializeHyperPlane(link_independent_generators);
                // O.initializeHyperPlane(link_independent_generators);
            }
            catch (int errorCode) {
                WARNING_PRINT("        CUDA & C++: Error initializing collision checking hyperplanes! Check previous error message!");
                throw;
            }

            auto stop1 = std::chrono::high_resolution_clock::now();
            auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(stop1 - start1);
            cout << "        CUDA & C++: Time taken by generating reachable sets: " << duration1.count() << " milliseconds" << endl;

            auto start2 = std::chrono::high_resolution_clock::now();

            // Solve optimization
            SmartPtr<armtd_NLP> mynlp = new armtd_NLP();
            try {
                mynlp->set_parameters(q_des, t_plan, &traj, &kd, &torque_radius, O_ptr.get(), u_s, surf_rad);
            }
            catch (int errorCode) {
                WARNING_PRINT("        CUDA & C++: Error initializing Ipopt! Check previous error message!");
                throw;
            }

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
            auto outputstream1 = std::ofstream(data_dir + ofile1);
            if (!outputstream1)
            {
                throw std::runtime_error("Open filestream failed!");
            }
            std::cout << "is_open: " << (bool)outputstream1 << std::endl;
            status = app->Initialize();
            if( status != Solve_Succeeded ) {
                WARNING_PRINT("Error during initialization!");
                outputstream1 << -1 << '\n';
                outputstream1.close();
                throw;
            }

            try {
                // Ask Ipopt to solve the problem
                status = app->OptimizeTNLP(mynlp);
            }
            catch (int errorCode) {
                WARNING_PRINT("        CUDA & C++: Error solving optimization problem! Check previous error message!");
                throw;
            }
            
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

            // set precision to 10 decimal digits
            outputstream1 << std::setprecision(10);

            // output k_opt
            auto k_opt = py::array_t<double>(NUM_FACTORS);
            double *k_opt_ptr = static_cast<double *>(k_opt.request().ptr);

            cout << "cpp planner raw solution: \n [ ";
            if (mynlp->feasible) {
                for (int i = 0; i < NUM_FACTORS; i++) {
                    cout << mynlp->solution[i] << ' ';
                    k_opt_ptr[i] = mynlp->solution[i] * k_range[i];
                }
            }
            else {
                for (int i = 0; i < NUM_FACTORS; i++) {
                    k_opt_ptr[i] = 0;
                }
                cout << "No feasible solution. ";
            }
            cout << "]" << endl;

            // output time cost (in milliseconds) in C++
            outputstream1 << duration1.count() + duration2.count();
            outputstream1.close();

            // output FRS and other information, you can comment them if they are unnecessary
            cout << "saving reach sets ..." << endl;
            std::ofstream outputstream2(data_dir + ofile2);
            outputstream2 << std::setprecision(10);
            for (int i = 0; i < NUM_TIME_STEPS; i++) {
                for (int j = 0; j < NUM_JOINTS; j++) {
                    for (int l = 0; l < 3; l++) {
                        outputstream2 << mynlp->link_sliced_center[i * NUM_JOINTS + j](l) << ' ';
                    }
                    outputstream2 << '\n';
                }
                outputstream2 << '\n';
            }
            outputstream2.close();

            std::ofstream outputstream3(data_dir + ofile3);
            outputstream3 << std::setprecision(10);
            for (int i = 0; i < NUM_TIME_STEPS; i++) {
                for (int j = 0; j < NUM_JOINTS; j++) {
                    for (int k = 0; k < 3; k++) {
                        for (int l = 0; l < 3 + 3; l++) {
                            outputstream3 << link_independent_generators[i * NUM_JOINTS + j](k, l) << ' ';
                        }
                        outputstream3 << '\n';
                    }
                    outputstream3 << '\n';
                }
                outputstream3 << '\n';
            }
            outputstream3.close();

            // std::ofstream outputstream4(outputfilename4);
            // outputstream4 << std::setprecision(10);
            // for (int i = 0; i < NUM_TIME_STEPS; i++) {
            //     for (int j = 0; j < NUM_FACTORS; j++) {
            //         outputstream4 << torque_radius(j, i) << ' '; // this is radius of final control input
            //     }
            //     outputstream4 << '\n';
            // }
            // outputstream4.close();

            // std::ofstream outputstream5(outputfilename5);
            // outputstream5 << std::setprecision(6);
            // for (int i = 0; i < mynlp->constraint_number; i++) {
            //     outputstream5 << mynlp->g_copy[i] << '\n';
            // }
            // outputstream5.close();
            
            return k_opt;
        }

        py::array_t<double> getDesTraj(Eigen::Ref<Eigen::VectorXd> q0, Eigen::Ref<Eigen::VectorXd> qd0, Eigen::Ref<Eigen::VectorXd> qdd0, py::array_t<double> k, double t){
            // return desired trajectory: [q, qd, qdd]
            int size = NUM_FACTORS * 3;
            auto traj_d = py::array_t<double>(size);
            double *traj_d_ptr = static_cast<double *>(traj_d.request().ptr);

            auto q0_ = q0;//.unchecked<1>();
            auto qd0_ = qd0 * DURATION;//.unchecked<1>();
            auto qdd0_ = qdd0 * DURATION * DURATION;//.unchecked<1>();
            t /= DURATION;
            auto k_ = k.unchecked<1>();

            for (int i = 0; i < NUM_FACTORS; i++){
                traj_d_ptr[3*i] = q_des_func(q0_(i), qd0_(i), qdd0_(i), k_(i), t);
                traj_d_ptr[3*i+1] = (qd_des_func(q0_(i), qd0_(i), qdd0_(i), k_(i), t))/DURATION;
                traj_d_ptr[3*i+2] = (qdd_des_func(q0_(i), qd0_(i), qdd0_(i), k_(i), t))/(DURATION*DURATION);
            }
            
            return traj_d;
        }

        // add collision checking function here and below in pybind function
        // take in the obstacles (from the pybind class above?), the precomputed joint locations from the matlab script, and adjacency matrix. 
        // return new adjacency matrix?

        

        // the pybind is initialized in bulletPlanner.py with the initialization of a armour planner

         py::array_t<double> graphCollisionCheck(Eigen::MatrixXd old_adj_matrix, Eigen::MatrixXd joint_positions) {
            /*
            Eigen::Ref<Eigen::MatrixXd>

            This function collision checks nodes in a graph and removes those in collision from the adjacency matrix. A new, collision-free adjacency
            matrix is returned.

            TODO: 
            - replace reading the old adjacency matrix from a file to taking in the adjacency matrix from python. can pass from csr_matrix?
            - can potentially return directly to csr_matrix? and also pass in from csr_matrix?
            */

            std::string output_message = "In C++ graphCollisionCheck";
            std::cout << output_message << std::endl;

            std::string pathname_collision = "/home/baiyuew/ROAHM/planning_wksp/src/kinova_perception/kinova_planning/src/rtd-pybullet/armour-dev/kinova_src/kinova_simulator_interfaces/kinova_planner_realtime/";

            // file names
            // composite set
                // joint_positions_composite.csv
                // adj_matrix_composite_range0p3.txt
                // edges: 28338492
            // uniform nodes hardware only version 2
                // joint_positions_uniform_v2.csv
                // adj_matrix_uniform_hardware_only_v2_range0p3.txt
                // edges: 16140775
            // uniform nodes hardware only version 3
                // joint_positions_uniform_v3.csv
                // adj_matrix_uniform_hardware_only_v3_range0p3.txt
                // edges: 17740918
            // uniform nodes hardware only version 4
                // joint_positions_uniform_v4.csv
                // adj_matrix_uniform_hardware_only_v4_range0p3.txt
                // edges: 17707982
            // uniform nodes hardware only version 4 dense and shifted valid only
                // joint_positions_uniform_dense_validOnly_v4.csv
                // adj_matrix_uniform_hardware_dense_shiftedValidOnly_v4_range0p3.txt
                // edges: 5482178
            // uniform nodes hardware dense shifted and random nodes
                // joint_positions_uniform_hardware_dense_rand.csv
                // adj_matrix_uniform_hardware_dense_rand_range0p3.txt
                // edges: 11930750
            // matlab files
                // joint_positions_uniform.csv
                // adj_matrix_uniform_mult5.txt

            // Hard Coded Files (need to be replaced if graph changes)
            const std::string inputfilename2 = pathname_collision + "joint_positions_uniform.csv";
            const std::string inputfilename3 = pathname_collision + "adj_matrix_uniform_mult5.csv";
            // const std::string outputfilename1 = pathname_collision + "node_feasibility.csv";
            // const std::string outputfilename2 = pathname_collision + "link_c.csv";
            const std::string outputfilename3 = pathname_collision + "collision_free_adj_matrix.csv";

            // O.initialize(obstacles, num_obstacles); // done above

            Eigen::Vector3d link_sliced_center[NUM_NODES_AT_ONE_TIME * NUM_JOINTS];
            Eigen::Matrix<double, 3, LINK_FRS_GENERATOR_NUM> link_independent_generators[NUM_NODES_AT_ONE_TIME * NUM_JOINTS];
            double* link_c = new double[NUM_JOINTS * NUM_NODES_AT_ONE_TIME * num_obstacles];
            bool* node_feasibilities = new bool[NUM_NODES];
            // std::ifstream inputstream2(inputfilename2);
            // std::ofstream outputstream1(outputfilename1);
            // std::ofstream outputstream2(outputfilename2);


            // pre-allocate new adjacency matrix

            // Eigen::MatrixXi new_adj_nodes = Eigen::MatrixXi::Zero(NUM_EDGES,2);
            // Eigen::VectorXd new_edge_distance = Eigen::VectorXd::Zero(NUM_EDGES);

            // Eigen::MatrixXd new_adj_matrix = Eigen::MatrixXd::Zero(NUM_EDGES,3);

            auto new_adj_matrix = py::array_t<double>(NUM_EDGES*3);
            double *new_adj_matrix_ptr = static_cast<double *>(new_adj_matrix.request().ptr);


            auto start1 = std::chrono::high_resolution_clock::now();

            int joint_positions_node_spacing = 9;
            int jpns = joint_positions_node_spacing;
            // reading joint positions
            for (int k = 0; k < NUM_NODES / NUM_NODES_AT_ONE_TIME; k++) {
                // int offset1 = k*NUM_NODES_AT_ONE_TIME-2*k;
                int set_offset = k*NUM_NODES_AT_ONE_TIME;
                for (int i = 0; i < NUM_NODES_AT_ONE_TIME; i++) {
                    Eigen::Vector3d pos1;
                    // inputstream2 >> pos1(0) >> pos1(1) >> pos1(2);
                    // extract the base joint position (every 8th in joint positions?)
                    // pos1(0) = joint_positions(i*jpns+offset1,0);
                    // pos1(1) = joint_positions(i*jpns+offset1,1);
                    // pos1(2) = joint_positions(i*jpns+offset1,2);
                    int first_row = (i+set_offset)*jpns;
                    pos1(0) = joint_positions(first_row,0);
                    pos1(1) = joint_positions(first_row,1);
                    pos1(2) = joint_positions(first_row,2);
                    Eigen::Vector3d pos2;
                    for (int j = 0; j < NUM_JOINTS; j++) {
                        // inputstream2 >> pos2(0) >> pos2(1) >> pos2(2);
                        // pos2(0) = joint_positions(i*jpns+j+offset1,0);
                        // pos2(1) = joint_positions(i*jpns+j+offset1,1);
                        // pos2(2) = joint_positions(i*jpns+j+offset1,2);
                        int curr_row = first_row+j+1;
                        pos2(0) = joint_positions(curr_row,0);
                        pos2(1) = joint_positions(curr_row,1);
                        pos2(2) = joint_positions(curr_row,2);

                        // if (k == 0){
                        //     if (i == 0){
                        //         std::cout << pos1(0) << ' ' << pos1(1) << ' ' << pos1(2) << std::endl;
                        //         std::cout << pos2(0) << ' ' << pos2(1) << ' ' << pos2(2) << std::endl;
                        //         std::cout << first_row << ' ' << curr_row << std::endl;
                        //         std::cout << std::endl;
                        //     }
                        // }
                        // if (k == 1){
                        //     if (i == 0){
                        //         std::cout << pos1(0) << ' ' << pos1(1) << ' ' << pos1(2) << std::endl;
                        //         std::cout << pos2(0) << ' ' << pos2(1) << ' ' << pos2(2) << std::endl;
                        //         std::cout << first_row << ' ' << curr_row << std::endl;
                        //         std::cout << std::endl;
                        //     }
                        // }
                        // if (k == 2){
                        //     if (i == 0){
                        //         std::cout << pos1(0) << ' ' << pos1(1) << ' ' << pos1(2) << std::endl;
                        //         std::cout << pos2(0) << ' ' << pos2(1) << ' ' << pos2(2) << std::endl;
                        //         std::cout << first_row << ' ' << curr_row << std::endl;
                        //         std::cout << std::endl;
                        //     }
                        // }
                        
                        // form the link zonotope
                        link_sliced_center[i * NUM_JOINTS + j] = 0.5 * (pos1 + pos2);
                        link_independent_generators[i * NUM_JOINTS + j] = 0.5 * (pos1 - pos2);
                        // TODO : ask Bohao why these indices are different than below loop
                        
                        // update previous position
                        pos1 = pos2;
                    }
                }

                /*
                Section II: Buffer obstacles and initialize collision checking hyperplanes
                */
                try {
                    O_ptr2->initializeHyperPlane(link_independent_generators);
                }
                catch (int errorCode) {
                    WARNING_PRINT("        CUDA & C++: Error initializing collision checking hyperplanes! Check previous error message!");
                    return new_adj_matrix; // -1; // 
                }

                /*
                Section III:
                    Collision checking
                */

                try {
                    O_ptr2->linkFRSConstraints(link_sliced_center, link_c);
                }
                catch (int errorCode) {std::
                    WARNING_PRINT("        CUDA & C++: Error peforming collision checking! Check previous error message!");
                    return new_adj_matrix; // -1; // 
                }

                /*
                Section IV:
                    Prepare output
                */
                // for (int i = 0; i < NUM_NODES_AT_ONE_TIME * num_obstacles; i++) {
                //     for (int j = 0; j < NUM_JOINTS; j++) {
                //         outputstream2 << link_c[j * NUM_NODES_AT_ONE_TIME * num_obstacles + i] << ' ';
                //     }
                //     outputstream2 << '\n';
                // }
                for (int i = 0; i < NUM_NODES_AT_ONE_TIME; i++) {
                    bool node_feasibility = true;
                    for (int j = 0; j < NUM_JOINTS; j++) {
                        for (int h = 0; h < num_obstacles; h++) {
                            if (link_c[(j * NUM_NODES_AT_ONE_TIME + i) * num_obstacles + h] > COLLISION_THRESHOLD) {
                                node_feasibility = false;
                                break;
                            }
                        }
                    }
                    // outputstream1 << node_feasibility << endl;
                    node_feasibilities[k * NUM_NODES_AT_ONE_TIME + i] = node_feasibility;
                }
            }

            auto stop1 = std::chrono::high_resolution_clock::now();
            auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(stop1 - start1);
            cout << "Time taken by peforming collision checking: " << duration1.count() << " milliseconds" << endl;

            // building adjacency matrix
            // std::ifstream inputstream3(inputfilename3);
            // std::ofstream outputstream3(outputfilename3);

            // for checking if specific nodes are feasible
            // if (!node_feasibilities[0]) {
            //     // checking 1st node for feasibility as start in hardware
            //     std::string output1 = "Start Node is infeasible";
            //     std::cout << std::endl << output1 << std::endl;
            // }
            // if (!node_feasibilities[1]) {
            //     // checking 2nd node for feasibility as goal in hardware
            //     std::string output2 = "Goal Node is infeasible";
            //     std::cout << std::endl << output2 << std::endl;
            // }
            // if (!node_feasibilities[13]) {
            //     // checking 20th node for feasibility as goal in hardware
            //     std::string output3 = "New Goal Node is infeasible";
            //     std::cout << std::endl << output3 << std::endl;
            // }
            
            auto start2 = std::chrono::high_resolution_clock::now();

            // build the new adjacency matrix
            for (int i = 0; i < NUM_EDGES; i++) {
                int node_a = old_adj_matrix(i,0);
                int node_b = old_adj_matrix(i,1);
                double edge_distance = old_adj_matrix(i,2);
                // inputstream3 >> node_a >> node_b >> edge_distance; // TODO replace this when passing in old adjacency matrix

                if (node_feasibilities[node_a] && node_feasibilities[node_b] && edge_distance < EDGE_THRESHOLD) {
                    // new_adj_nodes[i,0] = node_a;
                    // new_adj_nodes[i,1] = node_b;
                    // new_edge_distance[i,0] = edge_distance;

                    // std::cout << node_a << ' ' << node_b << ' ' << edge_distance << std::endl;
                    // std::cout << static_cast<float>(node_a) << ' ' << static_cast<float>(node_b) << ' ' << edge_distance << std::endl << std::endl;

                    // new_adj_matrix(i,0) = static_cast<float>(node_a);
                    // new_adj_matrix(i,1) = static_cast<float>(node_b); // trying to cast
                    // new_adj_matrix(i,2) = edge_distance;

                    new_adj_matrix_ptr[3*i] = float(node_a);
                    new_adj_matrix_ptr[3*i+1] = float(node_b); // trying to cast // _ptr
                    new_adj_matrix_ptr[3*i+2] = edge_distance;

                    // std::string output2 = "Feasible Node";
                    // std::cout << output2 << std::endl;

                    // outputstream3 << node_a << ' ' << node_b << ' ' << edge_distance << '\n';
                }
            }
            // could potentially parallelize the building of the adjacency matrix if not outputting to file?

            auto stop2 = std::chrono::high_resolution_clock::now();
            auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(stop2 - start2);
            cout << "Time taken by building new adjacency matrix: " << duration2.count() << " milliseconds" << endl;

            // inputstream2.close();
            // inputstream3.close();
            // outputstream1.close();
            // outputstream2.close();
            // outputstream3.close();
            delete[] link_c;
            delete[] node_feasibilities;

            return new_adj_matrix; // 0; // 

        }

        void free()
        {
            O_ptr.reset();
        }

        void set_datadir(const std::string &dir)
        {
            data_dir = dir;
        }

        std::string get_datadir()
        {
            return data_dir;
        }

};

PYBIND11_MODULE(armour_main_pybind, m) {
    py::class_<pzsparse>(m, "pzsparse")
        .def(py::init<py::array_t<double> &, const std::string &>())
        .def("getNumObstacles", &pzsparse::getNumObstacles)
        .def("optimize", &pzsparse::optimize)
        .def("getDesTraj", &pzsparse::getDesTraj)
        .def("setDataDir", &pzsparse::set_datadir)
        .def("getDataDir", &pzsparse::get_datadir)
        .def("graphCollisionCheck", &pzsparse::graphCollisionCheck, py::return_value_policy::copy); // reference_internal
        // .def("free", &pzsparse::free);
}