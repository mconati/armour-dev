#ifndef COLLISION_CHECKING_H
#define COLLISION_CHECKING_H

#include "Trajectory.h"

#define BUFFER_OBSTACLE_GENERATOR_NUM (MAX_OBSTACLE_GENERATOR_NUM + 3)
#define COMB_NUM BUFFER_OBSTACLE_GENERATOR_NUM * (BUFFER_OBSTACLE_GENERATOR_NUM - 1) / 2

__constant__ TYPE dev_obstacles[MAX_OBSTACLE_NUM * (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3];

__constant__ unsigned int dev_combA[COMB_NUM];
__constant__ unsigned int dev_combB[COMB_NUM];

class Obstacles {
public:
    int num_obstacles = 0;
    const TYPE* obstacles = nullptr;

    TYPE* dev_A = nullptr;
    TYPE* dev_d = nullptr;
    TYPE* dev_delta = nullptr;

	TYPE* dev_buffered_c = nullptr;
	TYPE* dev_buffered_G = nullptr;

    TYPE* dev_jointPositionRadius = nullptr;

    TYPE* dev_checkJointsPosition = nullptr;
    TYPE* dev_dk_checkJointsPosition = nullptr;

    TYPE* dev_link_c = nullptr;
    TYPE* dev_grad_link_c = nullptr;

	Obstacles(const TYPE* obstacles_inp, const int num_obstacles_inp);

	~Obstacles();

	void initializeHyperPlane(const TYPE* jointPositionRadius);

	void linkFRSConstraints(TYPE* checkJointsPosition,
                            TYPE* dk_checkJointsPosition,
							TYPE* link_c,
							TYPE* grad_link_c = nullptr);
};

__global__ void bufferObstaclesKernel(TYPE* jointPositionRadius, TYPE link_radius_x, TYPE link_radius_y, TYPE link_radius_z, TYPE* buffered_c, TYPE* buffered_G, int link_id);

__global__ void polytope_PH(TYPE* buffered_c, TYPE* buffered_G, TYPE* A, TYPE* d, TYPE* delta, unsigned int link_id);

__global__ void checkCollisionKernel(TYPE* A, TYPE* d, TYPE* delta, 
									 TYPE* checkJointsPosition,
									 TYPE* dk_checkJointsPosition,
									 int link_id, 
                                     TYPE* link_c, TYPE* grad_link_c);

#endif