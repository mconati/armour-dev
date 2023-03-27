#ifndef FAST_COLLISION_CHECKING_H
#define FAST_COLLISION_CHECKING_H

#include "Parameters.h"

#define NUM_NODES 6401000 // composite set of nodes: 6401000
#define NUM_NODES_AT_ONE_TIME 20000 // can not be larger than 20000 !!!
#define LINK_FRS_GENERATOR_NUM 1
#define BUFFER_OBSTACLE_GENERATOR_NUM_FAST (MAX_OBSTACLE_GENERATOR_NUM + LINK_FRS_GENERATOR_NUM)
#define COMB_NUM BUFFER_OBSTACLE_GENERATOR_NUM_FAST * (BUFFER_OBSTACLE_GENERATOR_NUM_FAST - 1) / 2

__constant__ double dev_obstacles_fast[MAX_OBSTACLE_NUM * (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3];

__constant__ uint dev_combA_fast[COMB_NUM];
__constant__ uint dev_combB_fast[COMB_NUM];

class SimplifiedObstacles {
public:
    int num_obstacles = 0;
    const double* obstacles = nullptr;

    Eigen::Vector3d* dev_A = nullptr;
    double* d = nullptr;
    double* dev_d = nullptr;
    double* delta = nullptr;
    double* dev_delta = nullptr;

	Eigen::Vector3d* dev_buffered_c = nullptr;
	Eigen::Matrix<double, 3, BUFFER_OBSTACLE_GENERATOR_NUM_FAST>* dev_buffered_G = nullptr;

    Eigen::Matrix<double, 3, LINK_FRS_GENERATOR_NUM>* dev_link_independent_generators = nullptr;

    Eigen::Vector3d* dev_link_sliced_center = nullptr;

    double* dev_link_c = nullptr;

    SimplifiedObstacles();

	~SimplifiedObstacles();

    void initialize(const double* obstacles_inp, const int num_obstacles_inp);

	void initializeHyperPlane(const Eigen::Matrix<double, 3, LINK_FRS_GENERATOR_NUM>* link_independent_generators);

	void linkFRSConstraints(Eigen::Vector3d* link_sliced_center,
							double* link_c);
};

__global__ void bufferObstaclesKernel(const Eigen::Matrix<double, 3, LINK_FRS_GENERATOR_NUM>* link_independent_generators, 
                                      Eigen::Vector3d* buffered_c, 
                                      Eigen::Matrix<double, 3, BUFFER_OBSTACLE_GENERATOR_NUM_FAST>* buffered_G, 
                                      int link_id);

__global__ void polytope_PH(Eigen::Vector3d* buffered_c, 
                            Eigen::Matrix<double, 3, BUFFER_OBSTACLE_GENERATOR_NUM_FAST>* buffered_G, 
                            Eigen::Vector3d* A, 
                            double* d, 
                            double* delta,
                            int link_id);

__global__ void checkCollisionKernel(Eigen::Vector3d* A, 
                                     double* d, 
                                     double* delta, 
									 Eigen::Vector3d* link_sliced_center,
									 int link_id, 
                                     double* link_c);

#endif
