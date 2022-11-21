#ifndef KINOVA_INFO_H
#define KINOVA_INFO_H

#include "Headers.h"

// Hard code robot physical properties here
// Kinova 7 dof

// number of joints
#define NUM_JOINTS 7

// number of factors (variables) in polynomial / number of actuated motors
// The last one joint (end effector) of Kinova is fixed
#define NUM_FACTORS 7

// 1,2,3 -> x,y,z, negative sign means rotate in reverse direction, 0 means fixed joint
const int axes[NUM_JOINTS] = { 3,3,3,3,3,3,3 }; 

// joint position translation element w.r.t previous joint frame, same as xyz in urdf
const TYPE trans[(NUM_JOINTS + 1) * 3] = {  0,          0,     0.15643,
											0,   0.005375,    -0.12838,
											0,   -0.21038,   -0.006375,
											0,   0.006375,    -0.21038,
											0,   -0.20843,   -0.006375,
											0, 0.00017505,    -0.10593,
											0,   -0.10593, -0.00017505,
											0,          0,           0 };


// joint position rotation element w.r.t previous joint frame, same as rpy in urdf
const TYPE rots[NUM_JOINTS * 3] = { M_PI,       0, 0,
                                    M_PI * 0.5, 0, 0,
								   -M_PI * 0.5, 0, 0,
								    M_PI * 0.5, 0, 0,
								   -M_PI * 0.5, 0, 0,
								    M_PI * 0.5, 0, 0,
								   -M_PI * 0.5, 0, 0 };

// link mass
const TYPE mass[NUM_JOINTS] = { 1.3773, 1.1636, 1.1636, 0.9302, 0.6781, 0.6781, 0.5 }; 
const TYPE mass_uncertainty = 0.03;

// link center of mass
const TYPE com[NUM_JOINTS * 3] = {  -0.000023, -0.010364,  -0.07336,
									-0.000044,  -0.09958, -0.013278,
									-0.000044, -0.006641, -0.117892,
									-0.000018, -0.075478, -0.015006,
									 0.000001, -0.009432, -0.063883,
									 0.000001, -0.045483,  -0.00965,
									 0.000281,  0.011402, -0.029798 };
const TYPE com_uncertainty = 0.0;

// link inertia
const TYPE inertia[NUM_JOINTS * 9] = {  0.00457, 0.000001, 0.000002, 0.000001, 0.004831, 0.000448, 0.000002, 0.000448, 0.001409,
										0.011088, 0.000005, 0, 0.000005, 0.001072, -0.000691, 0, -0.000691, 0.011255,
										0.010932, 0, -0.000007, 0, 0.011127, 0.000606, -0.000007, 0.000606, 0.001043,
										0.008147, -0.000001, 0, -0.000001, 0.000631, -0.0005, 0, -0.0005, 0.008316,
										0.001596, 0, 0, 0, 0.001607, 0.000256, 0, 0.000256, 0.000399,
										0.001641, 0, 0, 0, 0.00041, -0.000278, 0, -0.000278, 0.001641,
										0.000587, 0.000003, 0.000003, 0.000003, 0.000369, -0.000118, 0.000003, -0.000118, 0.000609 };
const TYPE inertia_uncertainty = 0.03;

// joint friction
// const TYPE friction[NUM_JOINTS] = {0.5217383101288284, 0.5769579059927288, 0.4213397946418778, 0.4945515376566732, 0.1611070502661354, 0.1333491185514130, 0.1434440181717370};
const TYPE friction[NUM_JOINTS] = {0.0}; // disable friction in Matlab simulation for now

// joint damping
// const TYPE damping[NUM_JOINTS] = {10.5, 7.4064845817230722, 9.9727633408172860, 8.2667950822503915, 8.8572249026528151, 8.7110831569332845, 8.8881903638306934};
const TYPE damping[NUM_JOINTS] = {0.0}; // disable damping in Matlab simulation for now

// joint armature / motor transmission inertia
const TYPE armature[NUM_JOINTS] = {8.03, 11.9962024615303644, 9.0025427861751517, 11.5806439316706360, 8.4665040917914123, 8.8537069373742430, 8.8587303664685315};

// other robot info
// 1000.0 means they are continuous joints and have no limits
const TYPE state_limits_lb[NUM_FACTORS] = { -1000.0,  -2.41,  -1000.0,  -2.66,  -1000.0,  -2.23,  -1000.0 }; // rad
const TYPE state_limits_ub[NUM_FACTORS] = {  1000.0,   2.41,   1000.0,   2.66,   1000.0,   2.23,   1000.0 }; // rad

const TYPE speed_limits[NUM_FACTORS] = { 1.3963, 1.3963, 1.3963, 1.3963, 1.2218, 1.2218, 1.2218 }; // rad/s

const TYPE torque_limits[NUM_FACTORS] = { 56.7, 56.7, 56.7, 56.7, 29.4, 29.4, 29.4 }; // N*m These are WARNING torque limits

const TYPE gravity = 9.81;

const TYPE link_radius[NUM_FACTORS][3] = {{0.0675, 0.0675, 0.0675}, 
									      {0.0675, 0.0675, 0.0675}, 
										  {0.0675, 0.0675, 0.0675}, 
										  {0.0675, 0.0675, 0.0675}, 
									      {0.0675, 0.0675, 0.0675}, 
										  {0.0675, 0.0675, 0.0675}, 
										  {0.10, 0.10, 0.10}};

// ultimate bound
const TYPE alpha = 1.0;
const TYPE V_m = 5e-5;
const TYPE M_max = 15.79635774;
const TYPE M_min = 5.095620491878957;
const TYPE eps = sqrt(2 * V_m / M_min);
const TYPE K = 10.0;
const TYPE qe = eps / K;
const TYPE qde = 2 * eps;
const TYPE qdae = eps;
const TYPE qddae = 2 * K * eps;

// We define this in order to deal with the fixed joints at the end (end effector) in kinova
// The number of 1 should be strictly equal to NUM_FACTORS !!!
const bool JOINTS_WE_CARE_IN_COLLISION_AVOIDANCE[NUM_JOINTS] = {1, 1, 1, 1, 1, 1, 1};

#endif
