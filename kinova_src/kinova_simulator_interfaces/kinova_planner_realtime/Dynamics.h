#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "Trajectory.h"

class KinematicsDynamics {
public:
	// orientation of frame i with respect to frame i-1
	matPZsparse R[NUM_JOINTS + 1];

	// orientation of frame i-1 with respect to frame i
    matPZsparse R_t[NUM_JOINTS];

	KinematicsDynamics() {};

	KinematicsDynamics(PZsparse* cosElt_arr, PZsparse* sinElt_arr);

	~KinematicsDynamics() {};

	void fk(PZsparse* p);

	void rnea(PZsparse* v_arr,
			  PZsparse* v_aux_arr,
			  PZsparse* a_arr,
			  PZsparse* mass_arr,
			  matPZsparse* I_arr,
			  PZsparse* u,
			  vecPZsparse f;
        	  vecPZsparse n;
			  bool setGravity = false);
};

#endif
