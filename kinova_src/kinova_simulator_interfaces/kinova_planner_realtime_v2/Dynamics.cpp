#ifndef DYNAMICS_CPP
#define DYNAMICS_CPP

#include "Dynamics.h"

KinematicsDynamics::KinematicsDynamics(const BezierCurve* traj_input) {
    traj = traj_input;

    links = PZsparseArray(NUM_FACTORS * 3, NUM_TIME_STEPS);
    mass_arr = PZsparseArray(NUM_JOINTS, 1);
    I_nominal_arr = PZsparseArray(NUM_JOINTS, 1);
    I_uncertain_arr = PZsparseArray(NUM_JOINTS, 1);
    u_nom = PZsparseArray(NUM_FACTORS, NUM_TIME_STEPS);
    u_nom_int = PZsparseArray(NUM_FACTORS, NUM_TIME_STEPS);
    r = PZsparseArray(NUM_FACTORS, 1);
    Mr = PZsparseArray(NUM_FACTORS, NUM_TIME_STEPS);
    // V = PZsparseArray(NUM_TIME_STEPS, 1);

    for (int i = 0; i < NUM_JOINTS; i++) {
        Eigen::MatrixXd mass_matrix(1, 1);
        mass_matrix(0) = mass[i];
        mass_arr(i) = PZsparse(mass_matrix, mass_uncertainty);

        Eigen::Matrix3d inertia_matrix;
        for (int j = 0; j < 9; j++) {
            inertia_matrix(j) = inertia[i * 9 + j]; // This may not be right...
        }
        I_nominal_arr(i) = PZsparse(inertia_matrix);
        I_uncertain_arr(i) = PZsparse(inertia_matrix, inertia_uncertainty);

        if (i < NUM_FACTORS) {
            r(i) = PZsparse(0, Interval(-eps, eps));
        }
    }
}

void KinematicsDynamics::fk() {
    PZsparse FK_R(3, 3);
    PZsparse FK_T(3, 1);
    int j = 0;

    for (int i = 0; i < NUM_JOINTS; i++) {
        Eigen::MatrixXd trans_matrix(3, 1);
        trans_matrix << trans[3 * i], trans[3 * i + 1], trans[3 * i + 2];
        PZsparse P(trans_matrix);
        
        FK_T = FK_T + FK_R * P;
        // FK_R = FK_R * R[i];
        
        // if (JOINTS_WE_CARE_IN_COLLISION_AVOIDANCE[i]) {
        //     p[j * 3 + 0] = *(FK_T.elt[0]);
        //     p[j * 3 + 1] = *(FK_T.elt[1]);
        //     p[j * 3 + 2] = *(FK_T.elt[2]);
        //     j++;
        // }
    }
}

// void KinematicsDynamics::rnea(PZsparseArray* R,
//                               PZsparseArray* R_t,
//                               PZsparseArray* v_arr,
//                               PZsparseArray* v_aux_arr,
//                               PZsparseArray* a_arr,
//                               PZsparseArray* mass_arr,
//                               PZsparseArray* I_arr,
//                               PZsparseArray* u,
//                               bool setGravity) {
//     PZsparse w(3, 1);
//     PZsparse wdot(3, 1);
//     PZsparse w_aux(3, 1);
//     PZsparse linear_acc(3, 1);

//     PZsparseArray F(NUM_JOINTS, 1);
//     PZsparseArray N(NUM_JOINTS, 1);

//     if (setGravity) { // set gravity
//         // directly modify the center of the PZ instance
//         linear_acc.center(2) = gravity;
//     }

//     // RNEA forward recursion
//     for (int i = 0; i < NUM_JOINTS; i++) {
//         // NOTE:
//         // This is just a simplified implementation!!!
//         // We assume all fixed joints are at the end and the revolute joints are consecutive
//         if (axes[i] != 0) { // revolute joints
//             // line 16
//             linear_acc = R_t[i] * (linear_acc 
//                                     + cross(wdot, trans + i * 3) 
//                                     + cross(w, cross(w_aux, trans + i * 3)));

//             // line 13
//             w = R_t[i] * w;
//             if (v_arr != nullptr) { // non-zero velocity input
//                 *(w.elt[abs(axes[i]) - 1]) += v_arr[i];
//             }

//             // line 14
//             w_aux = R_t[i] * w_aux;

//             // line 15
//             wdot = R_t[i] * wdot;

//             vecPZsparse temp; // temp = joint_vel(robot_params.q_index(i))*z(:,i)
//             if (v_arr != nullptr) { // non-zero velocity input
//                 *(temp.elt[abs(axes[i]) - 1]) = v_arr[i];
//             }

//             wdot = wdot + cross(w_aux, temp);

//             if (a_arr != nullptr) { // non-zero acceleration input
//                 *(wdot.elt[abs(axes[i]) - 1]) += a_arr[i];
//             }

//             // line 14
//             if (v_aux_arr != nullptr) { // non-zero auxiliary velocity input
//                 *(w_aux.elt[abs(axes[i]) - 1]) += v_aux_arr[i];
//             }
//         }
//         else { // fixed joints
//             // line 16
//             linear_acc = R_t[i] * (linear_acc 
//                                     + cross(wdot, trans + i * 3) 
//                                     + cross(w, cross(w_aux, trans + i * 3)));

//             // line 13
//             w = R_t[i] * w;

//             // line 14
//             w_aux = R_t[i] * w_aux;

//             // line 15
//             wdot = R_t[i] * wdot;
//         }

//         // line 23 & 27
//         if (mass_arr != nullptr) { // uncertain mass parameter
//             F[i] = mass_arr[i] * (linear_acc 
//                                     + cross(wdot, com + i * 3)
//                                     + cross(w, cross(w_aux, com + i * 3)));
//         }
//         else { // nominal mass parameter
//             F[i] = mass[i] * (linear_acc 
//                                 + cross(wdot, com + i * 3)
//                                 + cross(w, cross(w_aux, com + i * 3)));
//         }

//         // line 29
//         N[i] = I_arr[i] * wdot
//                 + cross(w_aux, (I_arr[i] * w));
//     }

//     PZsparse f(3, 1);
//     PZsparse n(3, 1);

//     // RNEA reverse recursion
//     for (int i = NUM_JOINTS - 1; i >= 0; i--) {
//         // line 29
//         n = N[i]
//             + R[i + 1] * n
//             + cross(com + i * 3, F[i])
//             + cross(trans + (i + 1) * 3, R[i + 1] * f);

//         // line 28
//         f = R[i + 1] * f + F[i];

//         if (axes[i] != 0) {
//             u[i] = *(n.elt[abs(axes[i]) - 1]);

//             u[i] = u[i] + armature[i] * a_arr[i];

//             u[i] = u[i] + damping[i] * v_arr[i];

//             // haven't implemented friction yet, this is more complicated...
//         }
//     }
// }

#endif