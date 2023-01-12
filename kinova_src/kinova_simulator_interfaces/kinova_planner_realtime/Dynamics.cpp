#ifndef DYNAMICS_CPP
#define DYNAMICS_CPP

#include "Dynamics.h"

KinematicsDynamics::KinematicsDynamics(PZsparse* cosElt_arr,
                                       PZsparse* sinElt_arr) {
    // Initialize rotation matrices
    for (int i = 0; i < NUM_JOINTS; i++) {
        if (sizeof(rots) == NUM_JOINTS * 3 * sizeof(TYPE)) { // non-trivial rpy information in joints
            R[i] = matPZsparse(rots[i * 3], rots[i * 3 + 1], rots[i * 3 + 2]);

            if (axes[i] != 0) { // revolute joints
                R[i] = R[i] * matPZsparse(cosElt_arr[i], sinElt_arr[i], axes[i]);
            }
        }
        else { // rpy all 0
            R[i] = matPZsparse(cosElt_arr[i], sinElt_arr[i], axes[i]);
        }

        transpose(R_t[i], R[i]);
    }

    // NOTE:
    // This is just a simplified implementation!!!
    // We always assume that number of bodies is equal to number of joints.
    // So We don't have to do anything with R[NUM_JOINTS], which is identity matrix by default
}

void KinematicsDynamics::fk(PZsparse* p) {
    matPZsparse FK_R;
    vecPZsparse FK_T;
    int j = 0;

    for (int i = 0; i < NUM_JOINTS; i++) {
        vecPZsparse P(trans + 3 * i);
        
        FK_T = FK_T + FK_R * P;
        FK_R = FK_R * R[i];
        
        if (JOINTS_WE_CARE_IN_COLLISION_AVOIDANCE[i]) {
            p[j * 3 + 0] = *(FK_T.elt[0]);
            p[j * 3 + 1] = *(FK_T.elt[1]);
            p[j * 3 + 2] = *(FK_T.elt[2]);
            j++;
        }
    }
}

void KinematicsDynamics::rnea(PZsparse* v_arr,
                              PZsparse* v_aux_arr,
                              PZsparse* a_arr,
                              PZsparse* mass_arr,
                              matPZsparse* I_arr,
                              PZsparse* u,
                              vecPZsparse* f_c,
                              vecPZsparse* n_c,
                              bool setGravity) {
    vecPZsparse w;
    vecPZsparse wdot;
    vecPZsparse w_aux;
    vecPZsparse linear_acc;

    vecPZsparse F[NUM_JOINTS];
    vecPZsparse N[NUM_JOINTS];

    if (setGravity) { // set gravity
        *(linear_acc.elt[2]) = PZsparse(gravity);
    }

    // RNEA forward recursion
    for (int i = 0; i < NUM_JOINTS; i++) {
        // NOTE:
        // This is just a simplified implementation!!!
        // We assume all fixed joints are at the end and the revolute joints are consecutive
        if (axes[i] != 0) { // revolute joints
            // line 16
            linear_acc = R_t[i] * (linear_acc 
                                    + cross(wdot, trans + i * 3) 
                                    + cross(w, cross(w_aux, trans + i * 3)));

            // line 13
            w = R_t[i] * w;
            if (v_arr != nullptr) { // non-zero velocity input
                *(w.elt[abs(axes[i]) - 1]) += v_arr[i];
            }

            // line 14
            w_aux = R_t[i] * w_aux;

            // line 15
            wdot = R_t[i] * wdot;

            vecPZsparse temp; // temp = joint_vel(robot_params.q_index(i))*z(:,i)
            if (v_arr != nullptr) { // non-zero velocity input
                *(temp.elt[abs(axes[i]) - 1]) = v_arr[i];
            }

            wdot = wdot + cross(w_aux, temp);

            if (a_arr != nullptr) { // non-zero acceleration input
                *(wdot.elt[abs(axes[i]) - 1]) += a_arr[i];
            }

            // line 14
            if (v_aux_arr != nullptr) { // non-zero auxiliary velocity input
                *(w_aux.elt[abs(axes[i]) - 1]) += v_aux_arr[i];
            }
        }
        else { // fixed joints
            // line 16
            linear_acc = R_t[i] * (linear_acc 
                                    + cross(wdot, trans + i * 3) 
                                    + cross(w, cross(w_aux, trans + i * 3)));

            // line 13
            w = R_t[i] * w;

            // line 14
            w_aux = R_t[i] * w_aux;

            // line 15
            wdot = R_t[i] * wdot;
        }

        // line 23 & 27
        if (mass_arr != nullptr) { // uncertain mass parameter
            F[i] = mass_arr[i] * (linear_acc 
                                    + cross(wdot, com + i * 3)
                                    + cross(w, cross(w_aux, com + i * 3)));
        }
        else { // nominal mass parameter
            F[i] = mass[i] * (linear_acc 
                                + cross(wdot, com + i * 3)
                                + cross(w, cross(w_aux, com + i * 3)));
        }

        // line 29
        N[i] = I_arr[i] * wdot
                + cross(w_aux, (I_arr[i] * w));
    }

    // pass out the contact joint through f_c and n_c
    vecPZsparse f;
    vecPZsparse n;

    // RNEA reverse recursion
    for (int i = NUM_JOINTS - 1; i >= 0; i--) {
        // line 29
        n = N[i]
            + R[i + 1] * n
            + cross(com + i * 3, F[i])
            + cross(trans + (i + 1) * 3, R[i + 1] * f);

        // line 28
        f = R[i + 1] * f + F[i];

        if (axes[i] != 0) {
            u[i] = *(n.elt[abs(axes[i]) - 1]);

            u[i] = u[i] + armature[i] * a_arr[i];

            u[i] = u[i] + damping[i] * v_arr[i];

            // haven't implemented friction yet, this is more complicated...
        }

        // passing out the contact joint RNEA outputs
        // ERROR! this needs to be NUM_JOINTS, i.e. 10th joint, but this for loop doesn't go that far.
        // I think the recursion is wrong? or the definition of NUM_JOINTS is wrong
        if (i == NUM_JOINTS-1 && f_c && n_c) {
            *f_c = f; // not sure how to assign these
            *n_c = n; // not sure how to assign these
        }
    }
}

#endif