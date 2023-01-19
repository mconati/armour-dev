#ifndef FETCH_INFO_CUH
#define FETCH_INFO_CUH

#include "Headers.h"

// Hard code robot physical properties here
// Fetch

#define NUM_JOINTS 9

// number of factors (variables) in polynomial
// The latter two joints of Fetch are fixed
#define NUM_FACTORS 7

// 1,2,3 -> x,y,z, negative sign means rotate in reverse direction, 0 means fixed joint
const int axes[NUM_JOINTS] = { 3,2,1,2,1,2,1,0,0 }; 

// joint position translation element w.r.t previous joint frame, same as xyz in urdf
const TYPE trans[(NUM_JOINTS + 1) * 3] = { -0.0326, 0, 0.726,
											 0.117, 0,  0.06,
											 0.219, 0,     0,
											 0.133, 0,     0,
										  	 0.197, 0,     0,
											0.1245, 0,     0,
											0.1385, 0,     0,
										   0.16645, 0,     0,
												 0, 0,     0,
												 0, 0,     0 };
												
// joint position rotation element w.r.t previous joint frame, same as rpy in urdf
// Fetch does not have a nontrivial rotation element for all joints!
const TYPE rots[1] = { 0 };

// link mass
const TYPE mass[NUM_JOINTS] = { 2.5587, 2.6615, 2.3311, 2.1299, 1.6563, 1.725, 0.1354, 1.5175, 2.26796 };
const TYPE mass_uncertainty = 0.03;

// link center of mass
const TYPE com[NUM_JOINTS * 3] = { 0.0927, -0.0056,  0.0564,
								   0.1432,  0.0072, -0.0001,
								   0.1165,  0.0014,       0,
								   0.1279,  0.0073,       0,
								   0.1097, -0.0266,       0,
								   0.0882,  0.0009, -0.0001,
								   0.0095,  0.0004, -0.0002,
								    -0.09, -0.0001, -0.0017,
								        0,       0,       0};
const TYPE com_uncertainty = 0.0;

// link inertia
const TYPE inertia[NUM_JOINTS * 9] = { 0.0043, -0.0001, 0.001, -0.0001, 0.0087, -0.0001, 0.001, -0.0001, 0.0087,
									   0.0028, -0.0021, 0, -0.0021, 0.0111, 0, 0, 0, 0.0112,
									   0.0019, -0.0001, 0, -0.0001, 0.0045, 0, 0, 0, 0.0047,
									   0.0024, -0.0016, 0, -0.0016, 0.0082, 0, 0, 0, 0.0084,
									   0.0016, -0.0003, 0, -0.0003, 0.003, 0, 0, 0, 0.0035,
									   0.0018, -0.0001, 0, -0.0001, 0.0042, 0, 0, 0, 0.0042,
									   0.0001, 0, 0, 0, 0.0001, 0, 0, 0, 0.0001,
									   0.0013, 0, 0, 0, 0.0019, 0, 0, 0, 0.0024,
									   0, 0, 0, 0, 0, 0, 0, 0, 0 };
const TYPE inertia_uncertainty = 0.03;

// joint friction
const TYPE friction[NUM_JOINTS] = {0};

// joint damping
const TYPE damping[NUM_JOINTS] = {0};

// joint armature / motor transmission inertia
const TYPE armature[NUM_JOINTS] = {0};

// ultimate bound
const TYPE V_m = 1e-7;
const TYPE M_min = 5.09562049;
const TYPE eps = sqrt(2 * V_m / M_min); // 0.0191;
const TYPE K = 5.0;
const TYPE qe = eps / K;
const TYPE qde = 2 * eps;
const TYPE qdae = eps;
const TYPE qddae = 2 * K * eps;

// robust controller info
const TYPE alpha = 1.0;

// other robot info
// 1000.0 means they are continuous joints and have no limits
const TYPE state_limits_lb[NUM_FACTORS] = { -1.6056,  -1.221,  -1000.0,  -2.251,  -1000.0,  -2.16,  -1000.0 };
const TYPE state_limits_ub[NUM_FACTORS] = {  1.6056,   1.518,   1000.0,   2.251,   1000.0,   2.16,   1000.0 };

const TYPE speed_limits[NUM_FACTORS] = { 1.256, 1.454, 1.571, 1.521, 1.571, 2.268, 2.268 };

const TYPE torque_limits[NUM_FACTORS] = { 33.82, 131.76, 76.94, 66.18, 29.35, 25.7, 7.36 };

const TYPE gravity = 9.81;

const TYPE link_radius[NUM_FACTORS][3] = {{0.04, 0.04, 0.04}, 
									      {0.04, 0.04, 0.04}, 
										  {0.04, 0.04, 0.04}, 
										  {0.04, 0.04, 0.04}, 
									      {0.04, 0.04, 0.04}, 
										  {0.04, 0.04, 0.04}, 
										  {0.04, 0.04, 0.04}};

// We define this in order to deal with the fixed joints at the end (end effector) in fetch
// The number of 1 should be strictly equal to NUM_FACTORS !!!
const bool JOINTS_WE_CARE_IN_COLLISION_AVOIDANCE[NUM_JOINTS] = {1, 1, 1, 1, 1, 1, 0, 1};

#endif
