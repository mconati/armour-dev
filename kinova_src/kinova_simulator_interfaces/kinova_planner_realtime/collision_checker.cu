#include "FastCollisionChecking.h"

const std::string inputfilename1 = "obstacles.csv";
const std::string inputfilename2 = "joint_positions.csv";
const std::string outputfilename1 = "node_feasibility.csv";
const std::string outputfilename2 = "link_c.csv";

#define COLLISION_THRESHOLD -0.05

int main() {
/*
Section I:
    Parse input
    There is no check and warning, so be careful!
*/

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

    SimplifiedObstacles O;
    Eigen::Vector3d link_sliced_center[NUM_NODES * NUM_JOINTS];
    Eigen::Matrix<double, 3, LINK_FRS_GENERATOR_NUM> link_independent_generators[NUM_NODES * NUM_JOINTS];

    const int num_obstacles = 10;
    double obstacles[MAX_OBSTACLE_NUM * (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3] = {0.0};

    std::ifstream inputstream1(inputfilename1);

    auto start0 = std::chrono::high_resolution_clock::now();
  
    // inputstream1 >> num_obstacles;
    if (num_obstacles > MAX_OBSTACLE_NUM || num_obstacles < 0) {
        throw;
    }
    if (num_obstacles > 0) {
        for (int i = 0; i < num_obstacles * (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3; i++) {
            inputstream1 >> obstacles[i];
        }
    }

    inputstream1.close();

    O.initialize(obstacles, num_obstacles); 

    std::ifstream inputstream2(inputfilename2);

    for (int i = 0; i < NUM_NODES; i++) {
        Eigen::Vector3d pos1;
        inputstream2 >> pos1(0) >> pos1(1) >> pos1(2);
        Eigen::Vector3d pos2;
        for (int j = 0; j < NUM_JOINTS; j++) {
            inputstream2 >> pos2(0) >> pos2(1) >> pos2(2);

            link_sliced_center[i * NUM_JOINTS + j] = 0.5 * (pos1 + pos2);
            link_independent_generators[i * NUM_JOINTS + j] = 0.5 * (pos1 - pos2);

            pos1 = pos2;
        }
    }

    inputstream2.close();

    auto stop0 = std::chrono::high_resolution_clock::now();
    auto duration0 = std::chrono::duration_cast<std::chrono::milliseconds>(stop0 - start0);
    cout << "        CUDA & C++: Time taken by loading data: " << duration0.count() << " milliseconds" << endl;

    double* link_c = new double[NUM_JOINTS * NUM_NODES * num_obstacles];
    std::ofstream outputstream1(outputfilename1);
    std::ofstream outputstream2(outputfilename2);

    for (int k = 0; k < NUM_NODES / NUM_NODES_AT_ONE_TIME; k++) {
        /*
        Section II: Buffer obstacles and initialize collision checking hyperplanes
        */
        auto start1 = std::chrono::high_resolution_clock::now();

        try {
            O.initializeHyperPlane(link_independent_generators + k * NUM_NODES_AT_ONE_TIME * NUM_JOINTS);
        }
        catch (int errorCode) {
            WARNING_PRINT("        CUDA & C++: Error initializing collision checking hyperplanes! Check previous error message!");
            return -1;
        }

        auto stop1 = std::chrono::high_resolution_clock::now();
        auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(stop1 - start1);
        // cout << "        CUDA & C++: Time taken by initializing collision checking hyperplanes: " << duration1.count() << " milliseconds" << endl;

        /*
        Section III:
            Collision checking
        */

        auto start2 = std::chrono::high_resolution_clock::now();

        try {
            O.linkFRSConstraints(link_sliced_center + k * NUM_NODES_AT_ONE_TIME * NUM_JOINTS, link_c);
        }
        catch (int errorCode) {
            WARNING_PRINT("        CUDA & C++: Error peforming collision checking! Check previous error message!");
            return -1;
        }
        
        auto stop2 = std::chrono::high_resolution_clock::now();
        auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(stop2 - start2);
        // cout << "        CUDA & C++: Time taken by peforming collision checking: " << duration1.count() << " milliseconds" << endl;

        /*
        Section IV:
            Prepare output
        */
        for (int i = 0; i < NUM_NODES_AT_ONE_TIME * num_obstacles; i++) {
            for (int j = 0; j < NUM_JOINTS; j++) {
                outputstream2 << link_c[j * NUM_NODES_AT_ONE_TIME * num_obstacles + i] << ' ';
            }
            outputstream2 << '\n';
        }
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
            outputstream1 << node_feasibility << endl;
        }
    }
    
    outputstream2.close();
    delete[] link_c;
    
    return 0;
}