#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "Trajectory.h"

class KinematicsDynamics {
public:
	const BezierCurve* traj = nullptr;

	// inertial paramete PZs
	// mass_nominal already stored as double array mass in RobotInfo.h
    PZsparseArray mass_arr;
    PZsparseArray I_nominal_arr;
    PZsparseArray I_uncertain_arr;

	// link PZs
	PZsparseArray links;

	// nominal torque PZs
    PZsparseArray u_nom;
    PZsparseArray u_nom_int;

	// other PZs
    PZsparseArray r;
    PZsparseArray Mr;
    // PZsparseArray V;

	KinematicsDynamics() {}

	KinematicsDynamics(const BezierCurve* traj_input);

	// generate link PZs through forward kinematics
	void fk();

	// generate nominal torque PZs through rnea
	// void rnea(PZsparseArray* R,
	// 		  PZsparseArray* R_t,
	// 		  PZsparseArray* v_arr,
	// 		  PZsparseArray* v_aux_arr,
	// 		  PZsparseArray* a_arr,
	// 		  PZsparseArray* mass_arr,
	// 		  PZsparseArray* I_arr,
	// 		  PZsparseArray* u,
	// 		  bool setGravity = true);
};

#endif
