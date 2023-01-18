#ifndef COLLISION_CHECKING_CU
#define COLLISION_CHECKING_CU

#include "CollisionChecking.h"

Obstacles::Obstacles(const TYPE* obstacles_inp, const int num_obstacles_inp) {
    obstacles = obstacles_inp;
    num_obstacles = num_obstacles_inp;

    if (num_obstacles > MAX_OBSTACLE_NUM) {
        WARNING_PRINT("Number of obstacles larger than MAX_OBSTACLE_NUM !\n");
        throw;
    }

    // build pre-defined hyper-plane equations for collision checking between links and obstacles
    if (obstacles != nullptr && num_obstacles > 0) {
        cudaMalloc((void**)&dev_jointPositionRadius, NUM_FACTORS * 3 * NUM_TIME_STEPS * sizeof(TYPE));
    
        cudaMalloc((void**)&dev_A, NUM_FACTORS * NUM_TIME_STEPS * num_obstacles * COMB_NUM * 3 * sizeof(TYPE));
        cudaMalloc((void**)&dev_d, NUM_FACTORS * NUM_TIME_STEPS * num_obstacles * COMB_NUM * sizeof(TYPE));
        cudaMalloc((void**)&dev_delta, NUM_FACTORS * NUM_TIME_STEPS * num_obstacles * COMB_NUM * sizeof(TYPE));

        // initialize constant memory
        // obstacle data
        cudaMemcpyToSymbol(dev_obstacles, obstacles, num_obstacles * (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3 * sizeof(TYPE));

        // combination index
        unsigned int combA[COMB_NUM], combB[COMB_NUM];
        unsigned int a_id = 0, b_id = 1;
        for (unsigned int i = 0; i < COMB_NUM; i++) {
            combA[i] = a_id;
            combB[i] = b_id;

            if (b_id < BUFFER_OBSTACLE_GENERATOR_NUM - 1) {
                b_id++;
            }
            else {
                a_id++;
                b_id = a_id + 1;
            }
        }

        cudaMemcpyToSymbol(dev_combA, combA, COMB_NUM * sizeof(unsigned int));
        cudaMemcpyToSymbol(dev_combB, combB, COMB_NUM * sizeof(unsigned int));

        cudaMalloc((void**)&dev_buffered_c, NUM_TIME_STEPS * num_obstacles * 3 * sizeof(TYPE));
        cudaMalloc((void**)&dev_buffered_G, NUM_TIME_STEPS * num_obstacles * (MAX_OBSTACLE_GENERATOR_NUM + 3) * 3 * sizeof(TYPE));

        cudaMalloc((void**)&dev_checkJointsPosition, NUM_TIME_STEPS * NUM_FACTORS * 3 * sizeof(TYPE));
        cudaMalloc((void**)&dev_dk_checkJointsPosition, NUM_TIME_STEPS * NUM_FACTORS * 3 * NUM_FACTORS * sizeof(TYPE));

        cudaMalloc((void**)&dev_link_c, NUM_TIME_STEPS * num_obstacles * sizeof(TYPE));
        cudaMalloc((void**)&dev_grad_link_c, NUM_TIME_STEPS * num_obstacles * NUM_FACTORS * sizeof(TYPE));
    }
}

Obstacles::~Obstacles() {
    cudaFree(dev_jointPositionRadius);

    cudaFree(dev_A);
    cudaFree(dev_d);
    cudaFree(dev_delta);

    cudaFree(dev_buffered_c);
    cudaFree(dev_buffered_G);

    cudaFree(dev_checkJointsPosition);
    cudaFree(dev_dk_checkJointsPosition);

    cudaFree(dev_link_c);
    cudaFree(dev_grad_link_c);
}

void Obstacles::initializeHyperPlane(const TYPE* jointPositionRadius) {
    if (num_obstacles == 0) return;

    cudaMemcpy(dev_jointPositionRadius, jointPositionRadius, NUM_TIME_STEPS * NUM_FACTORS * 3 * sizeof(TYPE), cudaMemcpyHostToDevice);

    for (unsigned int i = 0; i < NUM_FACTORS; i++) {
        cudaMemset(dev_buffered_G, 0, NUM_TIME_STEPS * num_obstacles * (MAX_OBSTACLE_GENERATOR_NUM + 3) * 3 * sizeof(TYPE));

        dim3 block1(3, num_obstacles);
        bufferObstaclesKernel << < NUM_TIME_STEPS, block1 >> > (dev_jointPositionRadius, link_radius[i][0], link_radius[i][1], link_radius[i][2], dev_buffered_c, dev_buffered_G, i);
        
        dim3 grid1(NUM_TIME_STEPS, num_obstacles);
        polytope_PH << < grid1, COMB_NUM >> > (dev_buffered_c, dev_buffered_G, dev_A, dev_d, dev_delta, i);
    }

    cudaDeviceSynchronize();
}

void Obstacles::linkFRSConstraints(TYPE* checkJointsPosition,
                                   TYPE* dk_checkJointsPosition,
                                   TYPE* link_c,
                                   TYPE* grad_link_c) {
    if (num_obstacles == 0) return;
                                    
    bool ifComputeGradient = true;

    cudaMemcpy(dev_checkJointsPosition, checkJointsPosition, NUM_TIME_STEPS * NUM_FACTORS * 3 * sizeof(TYPE), cudaMemcpyHostToDevice);

    if (dk_checkJointsPosition == nullptr || grad_link_c == nullptr) {
        ifComputeGradient = false;
    }
    else {
        cudaMemcpy(dev_dk_checkJointsPosition, dk_checkJointsPosition, NUM_TIME_STEPS * NUM_FACTORS * 3 * NUM_FACTORS * sizeof(TYPE), cudaMemcpyHostToDevice);
    }

	// ignore the first link (base -> first joint) collision checking
    for (unsigned int i = 1; i < NUM_FACTORS; i++) {
        dim3 grid1(NUM_TIME_STEPS, num_obstacles);
        
        if (ifComputeGradient) {
            checkCollisionKernel << < grid1, COMB_NUM >> > (dev_A, dev_d, dev_delta,
                                                            dev_checkJointsPosition,
                                                            dev_dk_checkJointsPosition,
                                                            i, 
                                                            dev_link_c, dev_grad_link_c);
        }
        else {
            checkCollisionKernel << < grid1, COMB_NUM >> > (dev_A, dev_d, dev_delta,
                                                            dev_checkJointsPosition,
                                                            nullptr,
                                                            i, 
                                                            dev_link_c, dev_grad_link_c);
        }
            
        if (link_c != nullptr) {
            cudaMemcpy(link_c + (i - 1) * NUM_TIME_STEPS * num_obstacles, dev_link_c, NUM_TIME_STEPS * num_obstacles * sizeof(TYPE), cudaMemcpyDeviceToHost);
        }
        if (grad_link_c != nullptr) {
            cudaMemcpy(grad_link_c + (i - 1) * NUM_TIME_STEPS * num_obstacles * NUM_FACTORS, dev_grad_link_c, NUM_TIME_STEPS * num_obstacles * NUM_FACTORS * sizeof(TYPE), cudaMemcpyDeviceToHost);
        }
    }
}

__global__ void bufferObstaclesKernel(TYPE* jointPositionRadius, TYPE link_radius_x, TYPE link_radius_y, TYPE link_radius_z, TYPE* buffered_c, TYPE* buffered_G, int link_id) {
	unsigned int time_id = blockIdx.x;
	unsigned int obs_id = threadIdx.y;
	unsigned int p_id = threadIdx.x;
	unsigned int num_obstacles = blockDim.y;

	__shared__ TYPE Grest[3];

	if (obs_id == 0) {
		if (p_id == 0) {
			Grest[0] = jointPositionRadius[(time_id * NUM_FACTORS + link_id) * 3    ] + link_radius_x;
		}
		else if (p_id == 1) {
			Grest[1] = jointPositionRadius[(time_id * NUM_FACTORS + link_id) * 3 + 1] + link_radius_y;
		}
		else {
			Grest[2] = jointPositionRadius[(time_id * NUM_FACTORS + link_id) * 3 + 2] + link_radius_z;
		}
	}

	__syncthreads();

	// copy obstacle center to buffered obstacle center
	buffered_c[(time_id * num_obstacles + obs_id) * 3 + p_id] = dev_obstacles[obs_id * (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3 + p_id];

	// copy obstacle generator to buffered obstacle generator
	for (unsigned int i = 0; i < MAX_OBSTACLE_GENERATOR_NUM; i++) {
		buffered_G[((time_id * num_obstacles + obs_id) * (MAX_OBSTACLE_GENERATOR_NUM + 3) + i) * 3 + p_id] = dev_obstacles[(obs_id * (MAX_OBSTACLE_GENERATOR_NUM + 1) + i + 1) * 3 + p_id];
	}

	// copy joint position error to buffered obstacle generator
	buffered_G[((time_id * num_obstacles + obs_id) * (MAX_OBSTACLE_GENERATOR_NUM + 3) + p_id + MAX_OBSTACLE_GENERATOR_NUM) * 3 + p_id] = Grest[p_id];
}

__global__ void polytope_PH(TYPE* buffered_c, TYPE* buffered_G, TYPE* A, TYPE* d, TYPE* delta, unsigned int link_id) {
	unsigned int time_id = blockIdx.x;
	unsigned int obs_id = blockIdx.y;
	unsigned int num_obstacles = gridDim.y;
	unsigned int p_id = threadIdx.x;

	unsigned int a_id = dev_combA[p_id];
	unsigned int b_id = dev_combB[p_id];

	__shared__ TYPE C[COMB_NUM][3];
	__shared__ TYPE G[3][MAX_OBSTACLE_GENERATOR_NUM + 3];
	__shared__ TYPE c[3];

	if (p_id < MAX_OBSTACLE_GENERATOR_NUM + 3) {
		G[0][p_id] = buffered_G[((time_id * num_obstacles + obs_id) * (MAX_OBSTACLE_GENERATOR_NUM + 3) + p_id) * 3    ];
		G[1][p_id] = buffered_G[((time_id * num_obstacles + obs_id) * (MAX_OBSTACLE_GENERATOR_NUM + 3) + p_id) * 3 + 1];
		G[2][p_id] = buffered_G[((time_id * num_obstacles + obs_id) * (MAX_OBSTACLE_GENERATOR_NUM + 3) + p_id) * 3 + 2];
	}

	if (p_id < 3) {
		c[p_id] = buffered_c[(time_id * num_obstacles + obs_id) * 3 + p_id];
	}

	__syncthreads();

	TYPE a_x = G[0][a_id];
	TYPE a_y = G[1][a_id];
	TYPE a_z = G[2][a_id];

	TYPE b_x = G[0][b_id];
	TYPE b_y = G[1][b_id];
	TYPE b_z = G[2][b_id];

	TYPE cross_x = a_y * b_z - a_z * b_y;
	TYPE cross_y = a_z * b_x - a_x * b_z;
	TYPE cross_z = a_x * b_y - a_y * b_x;

	TYPE norm_cross = sqrt(cross_x * cross_x + cross_y * cross_y + cross_z * cross_z);

	if (norm_cross > 0) {
		C[p_id][0] = cross_x / norm_cross;
		C[p_id][1] = cross_y / norm_cross;
		C[p_id][2] = cross_z / norm_cross;
	}
	else {
		C[p_id][0] = 0.0;
		C[p_id][1] = 0.0;
		C[p_id][2] = 0.0;
	}

	// A = C
	A[(((time_id * NUM_FACTORS + link_id) * num_obstacles + obs_id) * COMB_NUM + p_id) * 3    ] = C[p_id][0];
	A[(((time_id * NUM_FACTORS + link_id) * num_obstacles + obs_id) * COMB_NUM + p_id) * 3 + 1] = C[p_id][1];
	A[(((time_id * NUM_FACTORS + link_id) * num_obstacles + obs_id) * COMB_NUM + p_id) * 3 + 2] = C[p_id][2];

	// d = C * c
	d[((time_id * NUM_FACTORS + link_id) * num_obstacles + obs_id) * COMB_NUM + p_id] = C[p_id][0] * c[0] + C[p_id][1] * c[1] + C[p_id][2] * c[2];

	// delta = sum(abs(C * G))
	TYPE delta_res = 0.0;

	for (unsigned int j = 0; j < MAX_OBSTACLE_GENERATOR_NUM + 3; j++) {
		delta_res += fabs(C[p_id][0] * G[0][j] + C[p_id][1] * G[1][j] + C[p_id][2] * G[2][j]);
	}

	delta[((time_id * NUM_FACTORS + link_id) * num_obstacles + obs_id) * COMB_NUM + p_id] = delta_res;
}

__global__ void checkCollisionKernel(TYPE* A, TYPE* d, TYPE* delta, 
									 TYPE* checkJointsPosition,
									 TYPE* dk_checkJointsPosition,
									 int link_id, 
                                     TYPE* link_c, TYPE* grad_link_c) {
	unsigned int time_id = blockIdx.x;
	unsigned int obs_id = blockIdx.y;
	unsigned int num_obstacles = gridDim.y;
	unsigned int p_id = threadIdx.x;

	__shared__ TYPE lc1[3];
	__shared__ TYPE lc2[3];
	__shared__ TYPE d_lc1[3][NUM_FACTORS];
	__shared__ TYPE d_lc2[3][NUM_FACTORS];

	__shared__ TYPE C[COMB_NUM][3];

	__shared__ TYPE begin_dot[COMB_NUM];
	__shared__ TYPE end_dot[COMB_NUM];
	__shared__ TYPE delta_v[COMB_NUM];
	__shared__ TYPE lambda_lb[COMB_NUM];
	__shared__ TYPE lambda_ub[COMB_NUM];
	__shared__ TYPE LHS[COMB_NUM];

	if (p_id == 0) {
        if (link_id > 0) {
            lc1[0] = checkJointsPosition[(time_id * NUM_FACTORS + link_id - 1) * 3    ];
            lc1[1] = checkJointsPosition[(time_id * NUM_FACTORS + link_id - 1) * 3 + 1];
            lc1[2] = checkJointsPosition[(time_id * NUM_FACTORS + link_id - 1) * 3 + 2];
        }
        else {
            lc1[0] = 0; // base position x
            lc1[1] = 0; // base position y
            lc1[2] = 0; // base position z
        }
		lc2[0] = checkJointsPosition[(time_id * NUM_FACTORS + link_id) * 3    ];
		lc2[1] = checkJointsPosition[(time_id * NUM_FACTORS + link_id) * 3 + 1];
		lc2[2] = checkJointsPosition[(time_id * NUM_FACTORS + link_id) * 3 + 2];
	}
	
	if (p_id < NUM_FACTORS && dk_checkJointsPosition != nullptr) { // Assume that COMB_NUM is definitely larger than NUM_FACTORS
        if (link_id > 0) {
            d_lc1[0][p_id] = dk_checkJointsPosition[((time_id * NUM_FACTORS + link_id - 1) * 3    ) * NUM_FACTORS + p_id];
            d_lc1[1][p_id] = dk_checkJointsPosition[((time_id * NUM_FACTORS + link_id - 1) * 3 + 1) * NUM_FACTORS + p_id];
            d_lc1[2][p_id] = dk_checkJointsPosition[((time_id * NUM_FACTORS + link_id - 1) * 3 + 2) * NUM_FACTORS + p_id];
        }
        else {
            d_lc1[0][p_id] = 0;
            d_lc1[1][p_id] = 0;
            d_lc1[2][p_id] = 0;
        }
		d_lc2[0][p_id] = dk_checkJointsPosition[((time_id * NUM_FACTORS + link_id) * 3    ) * NUM_FACTORS + p_id];
		d_lc2[1][p_id] = dk_checkJointsPosition[((time_id * NUM_FACTORS + link_id) * 3 + 1) * NUM_FACTORS + p_id];
		d_lc2[2][p_id] = dk_checkJointsPosition[((time_id * NUM_FACTORS + link_id) * 3 + 2) * NUM_FACTORS + p_id];
	}

	lambda_lb[p_id] = -100000000;
	lambda_ub[p_id] = 100000000;
	LHS[p_id] = 0;

	__syncthreads();

	TYPE C1 = A[(((time_id * NUM_FACTORS + link_id) * num_obstacles + obs_id) * COMB_NUM + p_id) * 3    ];
	TYPE C2 = A[(((time_id * NUM_FACTORS + link_id) * num_obstacles + obs_id) * COMB_NUM + p_id) * 3 + 1];
	TYPE C3 = A[(((time_id * NUM_FACTORS + link_id) * num_obstacles + obs_id) * COMB_NUM + p_id) * 3 + 2];
	C[p_id][0] = C1;
	C[p_id][1] = C2;
	C[p_id][2] = C3;

	if (C1 != 0 || C2 != 0 || C3 != 0) {
		TYPE d_v = d[((time_id * NUM_FACTORS + link_id) * num_obstacles + obs_id) * COMB_NUM + p_id];
		delta_v[p_id] = delta[((time_id * NUM_FACTORS + link_id) * num_obstacles + obs_id) * COMB_NUM + p_id];
		begin_dot[p_id] = C1 * lc1[0] + C2 * lc1[1] + C3 * lc1[2] - d_v;
		end_dot[p_id] = C1 * lc2[0] + C2 * lc2[1] + C3 * lc2[2] - d_v;
		LHS[p_id] = begin_dot[p_id] - end_dot[p_id];

		if (LHS[p_id] > 0) {
			lambda_ub[p_id] = (delta_v[p_id] - end_dot[p_id]) / LHS[p_id];
            lambda_lb[p_id] = (-delta_v[p_id] - end_dot[p_id]) / LHS[p_id];
		}
		else {
			lambda_ub[p_id] = (-delta_v[p_id] - end_dot[p_id]) / LHS[p_id];
            lambda_lb[p_id] = (delta_v[p_id] - end_dot[p_id]) / LHS[p_id];
		}
	}

	__syncthreads();
	
	if (p_id == 0) {
		unsigned int lb_max_id = 0;
		unsigned int ub_min_id = 0;
		for (unsigned int i = 1; i < COMB_NUM; i++) {
			if (lambda_lb[lb_max_id] < lambda_lb[i]) {
				lb_max_id = i;
			}
			if (lambda_ub[i] < lambda_ub[ub_min_id]) {
				ub_min_id = i;
			}
		}
		
		const TYPE& lambda_lb_max = lambda_lb[lb_max_id];
		const TYPE& LHS_lb = LHS[lb_max_id];
		const TYPE begin_dot_lb = begin_dot[lb_max_id];
		const TYPE end_dot_lb = end_dot[lb_max_id];
		const TYPE& delta_lb = delta_v[lb_max_id];

		const TYPE& lambda_ub_min = lambda_ub[ub_min_id];
		const TYPE& LHS_ub = LHS[ub_min_id];
		// const TYPE begin_dot_ub = begin_dot[ub_min_id];
		const TYPE end_dot_ub = end_dot[ub_min_id];
		const TYPE& delta_ub = delta_v[ub_min_id];
		
		TYPE begin_dot_lb_derivative[NUM_FACTORS];
		TYPE end_dot_lb_derivative[NUM_FACTORS];
		TYPE begin_dot_ub_derivative[NUM_FACTORS];
		TYPE end_dot_ub_derivative[NUM_FACTORS];

		if (dk_checkJointsPosition != nullptr) {
			for (unsigned int i = 0; i < NUM_FACTORS; i++) {
				begin_dot_lb_derivative[i] = C[lb_max_id][0] * d_lc1[0][i] + C[lb_max_id][1] * d_lc1[1][i] + C[lb_max_id][2] * d_lc1[2][i];
				end_dot_lb_derivative[i] = C[lb_max_id][0] * d_lc2[0][i] + C[lb_max_id][1] * d_lc2[1][i] + C[lb_max_id][2] * d_lc2[2][i];
				begin_dot_ub_derivative[i] = C[ub_min_id][0] * d_lc1[0][i] + C[ub_min_id][1] * d_lc1[1][i] + C[ub_min_id][2] * d_lc1[2][i];
				end_dot_ub_derivative[i] = C[ub_min_id][0] * d_lc2[0][i] + C[ub_min_id][1] * d_lc2[1][i] + C[ub_min_id][2] * d_lc2[2][i];
			}
		}

		if (1 < lambda_lb_max) {
			if (LHS_lb > 0) {
				link_c[time_id * num_obstacles + obs_id] = begin_dot_lb - delta_lb;
				
				for (unsigned int i = 0; i < NUM_FACTORS; i++) {
					grad_link_c[(time_id * num_obstacles + obs_id) * NUM_FACTORS + i] = begin_dot_lb_derivative[i];
				}
			}
			else {
				link_c[time_id * num_obstacles + obs_id] = -begin_dot_lb - delta_lb;

				for (unsigned int i = 0; i < NUM_FACTORS; i++) {
					grad_link_c[(time_id * num_obstacles + obs_id) * NUM_FACTORS + i] = -begin_dot_lb_derivative[i];
				}
			}
		}
		else if (lambda_ub_min < 0) {
			if (LHS_ub > 0) {
				link_c[time_id * num_obstacles + obs_id] = -end_dot_ub + delta_ub;

				for (unsigned int i = 0; i < NUM_FACTORS; i++) {
					grad_link_c[(time_id * num_obstacles + obs_id) * NUM_FACTORS + i] = -end_dot_ub_derivative[i];
				}
			}
			else {
				link_c[time_id * num_obstacles + obs_id] = end_dot_ub + delta_ub;

				for (unsigned int i = 0; i < NUM_FACTORS; i++) {
					grad_link_c[(time_id * num_obstacles + obs_id) * NUM_FACTORS + i] = end_dot_ub_derivative[i];
				}
			}
		}
		else {
			TYPE LHS_ub_derivative[NUM_FACTORS];
			TYPE LHS_lb_derivative[NUM_FACTORS];

			if (dk_checkJointsPosition != nullptr) {
				for (unsigned int i = 0; i < NUM_FACTORS; i++) {
					LHS_ub_derivative[i] = begin_dot_ub_derivative[i] - end_dot_ub_derivative[i];
					LHS_lb_derivative[i] = begin_dot_lb_derivative[i] - end_dot_lb_derivative[i];
				}
			}

			if (LHS_lb > 0) {
				if (LHS_ub > 0) {
					link_c[time_id * num_obstacles + obs_id] = LHS_ub * (end_dot_lb + delta_lb) + LHS_lb * (-end_dot_ub + delta_ub);

					for (unsigned int i = 0; i < NUM_FACTORS; i++) {
						grad_link_c[(time_id * num_obstacles + obs_id) * NUM_FACTORS + i] = 
										(LHS_ub_derivative[i] * (end_dot_lb + delta_lb)) + 
										(LHS_ub * end_dot_lb_derivative[i]) + 
										(LHS_lb_derivative[i] * (-end_dot_ub + delta_ub)) + 
										(LHS_lb * (-end_dot_ub_derivative[i]));
					}
				}
				else {
					link_c[time_id * num_obstacles + obs_id] = -LHS_ub * (end_dot_lb + delta_lb) + LHS_lb * (end_dot_ub + delta_ub);

					for (unsigned int i = 0; i < NUM_FACTORS; i++) {
						grad_link_c[(time_id * num_obstacles + obs_id) * NUM_FACTORS + i] = 
										(-LHS_ub_derivative[i] * (end_dot_lb + delta_lb)) + 
										(-LHS_ub * end_dot_lb_derivative[i]) + 
										(LHS_lb_derivative[i] * (end_dot_ub + delta_ub)) + 
										(LHS_lb * (end_dot_ub_derivative[i]));
					}
				}
			}
			else {
				if (LHS_ub > 0) {
					link_c[time_id * num_obstacles + obs_id] = LHS_ub * (-end_dot_lb + delta_lb) - LHS_lb * (-end_dot_ub + delta_ub);

					for (unsigned int i = 0; i < NUM_FACTORS; i++) {
						grad_link_c[(time_id * num_obstacles + obs_id) * NUM_FACTORS + i] = 
										(LHS_ub_derivative[i] * (-end_dot_lb + delta_lb)) + 
										(LHS_ub * (-end_dot_lb_derivative[i])) + 
										(-LHS_lb_derivative[i] * (-end_dot_ub + delta_ub)) + 
										(LHS_lb * (end_dot_ub_derivative[i]));
					}
				}
				else {
					link_c[time_id * num_obstacles + obs_id] = -LHS_ub * (-end_dot_lb + delta_lb) - LHS_lb * (end_dot_ub + delta_ub);
				
					for (unsigned int i = 0; i < NUM_FACTORS; i++) {
						grad_link_c[(time_id * num_obstacles + obs_id) * NUM_FACTORS + i] = 
										(-LHS_ub_derivative[i] * (-end_dot_lb + delta_lb)) + 
										(LHS_ub * (end_dot_lb_derivative[i])) + 
										(-LHS_lb_derivative[i] * (end_dot_ub + delta_ub)) + 
										(-LHS_lb * (end_dot_ub_derivative[i]));
					}
				}
			}
		}
	}
}

#endif