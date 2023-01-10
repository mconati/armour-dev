#ifndef NLP_CLASS_CU
#define NLP_CLASS_CU

#include "NLPclass.h"

// constructor
armtd_NLP::armtd_NLP()
{
    checkJointsPosition = new TYPE[NUM_TIME_STEPS * NUM_FACTORS * 3];
    dk_checkJointsPosition = new TYPE[NUM_TIME_STEPS * NUM_FACTORS * 3 * NUM_FACTORS];
}


// destructor
armtd_NLP::~armtd_NLP()
{
    delete[] checkJointsPosition;
    delete[] dk_checkJointsPosition;
    delete[] g_copy;
}


bool armtd_NLP::set_parameters(
    TYPE* q_des_input,
    TYPE t_plan_input,
    BezierCurve* desired_trajectory_input,
    PZsparse* joint_position_input,
    PZsparse* control_input_input,
    TYPE* v_norm_input,
    Obstacles* obstacles_input,
    vecPZsparse* f_c_input,
    vecPZsparse* n_c_input,
    const Number* u_s,
    const Number* surf_rad
 ) 
 {
    q_des = q_des_input;
    t_plan = t_plan_input;
    desired_trajectory = desired_trajectory_input;
    joint_position = joint_position_input;
    control_input = control_input_input;
    v_norm = v_norm_input;
    obstacles = obstacles_input;
    f_c = f_c_input;
    n_c = n_c_input;

    constraint_number = NUM_FACTORS * NUM_TIME_STEPS +
                        (NUM_FACTORS - 1) * NUM_TIME_STEPS * obstacles->num_obstacles + 
                        NUM_FACTORS * 4 + 
                        NUM_TIME_STEPS * 3;
    // 1. torque constraints
    // 2. collision checking constraints (ignoring base link and no self-collision)
    // 3. joint limit constraints, pos/vel lower and upper constraints
    // 4. force constraints

    g_copy = new Number[constraint_number];

    return true;
}


bool armtd_NLP::get_nlp_info(
   Index&          n,
   Index&          m,
   Index&          nnz_jac_g,
   Index&          nnz_h_lag,
   IndexStyleEnum& index_style
)
{
    // The problem described NUM_FACTORS variables, x[NUM_FACTORS] through x[NUM_FACTORS] for each joint
    n = NUM_FACTORS;

    // number of inequality constraint
    m = constraint_number;

    nnz_jac_g = m * n;

    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}
// [TNLP_get_nlp_info]

// [TNLP_get_bounds_info]
// returns the variable bounds
bool armtd_NLP::get_bounds_info(
   Index   n,
   Number* x_l,
   Number* x_u,
   Index   m,
   Number* g_l,
   Number* g_u
)
{
    // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
    // If desired, we could assert to make sure they are what we think they are.
    if(n != NUM_FACTORS){
        WARNING_PRINT("*** Error wrong value of n in get_bounds_info!");
    }
    if(m != constraint_number){
        WARNING_PRINT("*** Error wrong value of m in get_bounds_info!");
    }

    // lower bounds
    for( Index i = 0; i < n; i++ ) {
        x_l[i] = -1.0;
    }

    // upper bounds  
    for( Index i = 0; i < n; i++ ) {
        x_u[i] = 1.0;
    }

    // control input constraints
    for( Index i = 0; i < NUM_TIME_STEPS; i++ ) {
        for( Index j = 0; j < NUM_FACTORS; j++ ) {
            g_l[i * NUM_FACTORS + j] = -torque_limits[j] + v_norm[i * NUM_FACTORS + j];
            g_u[i * NUM_FACTORS + j] = torque_limits[j] - v_norm[i * NUM_FACTORS + j];
        }
    }    
    Index offset = NUM_FACTORS * NUM_TIME_STEPS;

    // force constraints on contact joint between tray and cup
    // need to add the radius of the constraints here (like how is done with the torque constraints, v_norm passes out the radius I think)
    // double check signs and how the torque constraints are buffered using the radius
    // for v_norm in armour_main.cpp : getRadius(u_nom[t_ind * NUM_FACTORS + i].independent)

    //     separation constraint
    // upper bound should be zero and lower bound should be -inf
    for( Index i = offset; i < offset + NUM_TIME_STEPS; i++){
        g_l[i] = -1e19;
        g_u[i] = 0;
    }
    offset += NUM_TIME_STEPS;

    //     slipping constraint
    // upper bound should be zero and lower bound should be -inf for the reformulated constraint (not the normal friction law?)
    for( Index i = offset; i < offset + NUM_TIME_STEPS; i++){
        g_l[i] = -1e19;
        g_u[i] = 0;
    }
    offset += NUM_TIME_STEPS;

    //     tipping constraint
    // upper bound should be zero and lower bound should be -inf for the reformulated constraint
    for( Index i = offset; i < offset + NUM_TIME_STEPS; i++){
        g_l[i] = -1e19;
        g_u[i] = 0;
    }
    offset += NUM_TIME_STEPS;

    // collision avoidance constraints
    for( Index i = offset; i < offset + (NUM_FACTORS - 1) * NUM_TIME_STEPS * obstacles->num_obstacles; i++ ) {
        g_l[i] = -1e19;
        g_u[i] = 0;
    }
    offset += (NUM_FACTORS - 1) * NUM_TIME_STEPS * obstacles->num_obstacles;

    // state limit constraints
    //     minimum joint position
    for( Index i = offset; i < offset + NUM_FACTORS; i++ ) {
        g_l[i] = state_limits_lb[i - offset] + qe;
        g_u[i] = state_limits_ub[i - offset] - qe;
    }
    offset += NUM_FACTORS;

    //     maximum joint position
    for( Index i = offset; i < offset + NUM_FACTORS; i++ ) {
        g_l[i] = state_limits_lb[i - offset] + qe;
        g_u[i] = state_limits_ub[i - offset] - qe;
    }
    offset += NUM_FACTORS;

    //     minimum joint velocity
    for( Index i = offset; i < offset + NUM_FACTORS; i++ ) {
        g_l[i] = -speed_limits[i - offset] + qde;
        g_u[i] = speed_limits[i - offset] - qde;
    }
    offset += NUM_FACTORS;

    //     maximum joint velocity
    for( Index i = offset; i < offset + NUM_FACTORS; i++ ) {
        g_l[i] = -speed_limits[i - offset] + qde;
        g_u[i] = speed_limits[i - offset] - qde;
    }
    offset += NUM_FACTORS;

    return true;
}
// [TNLP_get_bounds_info]

// [TNLP_get_starting_point]
// returns the initial point for the problem
bool armtd_NLP::get_starting_point(
    Index   n,
    bool    init_x,
    Number* x,
    bool    init_z,
    Number* z_L,
    Number* z_U,
    Index   m,
    bool    init_lambda,
    Number* lambda
)
{
    // Here, we assume we only have starting values for x, if you code
    // your own NLP, you can provide starting values for the dual variables
    // if you wish
    if(init_x == false || init_z == true || init_lambda == true){
        WARNING_PRINT("*** Error wrong value of init in get_starting_point!");
    }

    if(n != NUM_FACTORS){
        WARNING_PRINT("*** Error wrong value of n in get_starting_point!");
    }

    for( Index i = 0; i < n; i++ ) {
        // initialize to zero
        // x[i] = 0.0;

        // try to avoid local minimum
        x[i] = min(max((q_des[i] - desired_trajectory->q0[i]) / k_range[i], -0.5), 0.5);
    }

    return true;
}
// [TNLP_get_starting_point]

// [TNLP_eval_f]
// returns the value of the objective function
bool armtd_NLP::eval_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number&       obj_value
)
{
    if(n != NUM_FACTORS){
       WARNING_PRINT("*** Error wrong value of n in eval_f!");
    }

    // obj_value = sum((q_plan - q_des).^2);
    obj_value = 0; 
    for(Index i = 0; i < n; i++){
        TYPE q_plan = q_des_func(desired_trajectory->q0[i], desired_trajectory->qd0[i], desired_trajectory->qdd0[i], k_range[i] * x[i], t_plan);
        obj_value += pow(q_plan - q_des[i], 2);
    }

    obj_value *= 100.0;

    return true;
}
// [TNLP_eval_f]

// [TNLP_eval_grad_f]
// return the gradient of the objective function grad_{x} f(x)
bool armtd_NLP::eval_grad_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number*       grad_f
)
{
    if(n != NUM_FACTORS){
        WARNING_PRINT("*** Error wrong value of n in eval_grad_f!");
    }

    for(Index i = 0; i < n; i++){
        TYPE q_plan = q_des_func(desired_trajectory->q0[i], desired_trajectory->qd0[i], desired_trajectory->qdd0[i], k_range[i] * x[i], t_plan);
        TYPE dk_q_plan = pow(t_plan,3) * (6 * pow(t_plan,2) - 15 * t_plan + 10);
        grad_f[i] = (2 * (q_plan - q_des[i]) * dk_q_plan * k_range[i]) * 100.0;
    }

    return true;
}
// [TNLP_eval_grad_f]

// [TNLP_eval_g]
// return the value of the constraints: g(x)
bool armtd_NLP::eval_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Number*       g
)
{
    if(n != NUM_FACTORS){
        WARNING_PRINT("*** Error wrong value of n in eval_g!");
    }
    if(m != constraint_number){
        WARNING_PRINT("*** Error wrong value of m in eval_g!");
    }

    Index i;
    #pragma omp parallel for private(i) schedule(static, (NUM_TIME_STEPS * NUM_FACTORS + NUM_TIME_STEPS) / NUM_THREADS)
    for(i = 0; i < (NUM_TIME_STEPS * NUM_FACTORS + NUM_TIME_STEPS); i++) {  // took out the times 3 here since I am calculating all three contact constraints at the same time, but they still fill up all the same space
        
        if(i < (NUM_TIME_STEPS * NUM_FACTORS)) {
            // Part 1. slice the control input PZ to get the center of the input bound, 
            //         while the radius of the input bound has already been incorporated in ipopt constraint bounds
            g[i] = getCenter(control_input[i].slice(x));

            // Part 2. slice the forward kinematics PZ to get the center of the joint position bound, 
            //         while the radius of the joint position bound has already been considered in Obstacles class
            checkJointsPosition[i * 3    ] = getCenter(joint_position[i * 3    ].slice(x));
            checkJointsPosition[i * 3 + 1] = getCenter(joint_position[i * 3 + 1].slice(x));
            checkJointsPosition[i * 3 + 2] = getCenter(joint_position[i * 3 + 2].slice(x));
        } 
        else {
            // Part 3. force constraints on contact joint between tray and object
            // note: could change the next couple of lines to not slice multiple times

            //     force
            // get centers and their squares
            Number f_c_x_center = getCenter(f_c[i].elt[0].slice(x));
            Number f_c_x_center_2 = f_c_x_center * f_c_x_center;
            Number f_c_y_center = getCenter(f_c[i].elt[1].slice(x));
            Number f_c_y_center_2 = f_c_y_center * f_c_y_center;
            Number f_c_z_center = getCenter(f_c[i].elt[2].slice(x));
            Number f_c_z_center_2 = f_c_z_center * f_c_z_center;
            // get radii of independent generators and their squares
            TYPE f_c_x_radius = getRadius(f_c[i].elt[0].slice(x).independent)
            TYPE f_c_x_radius_2 = f_c_x_radius * f_c_x_radius;
            TYPE f_c_y_radius = getRadius(f_c[i].elt[1].slice(x).independent)
            TYPE f_c_y_radius_2 = f_c_y_radius * f_c_y_radius;
            TYPE f_c_z_radius = getRadius(f_c[i].elt[2].slice(x).independent)
            TYPE f_c_z_radius_2 = f_c_z_radius * f_c_z_radius;

            //     moment
            // get centers and their squares
            Number n_c_x_center = getCenter(n_c[i].elt[0].slice(x));
            Number n_c_x_center_2 = n_c_x_center * n_c_x_center;
            Number n_c_y_center = getCenter(n_c[i].elt[1].slice(x));
            Number n_c_y_center_2 = n_c_y_center * n_c_y_center;
            Number n_c_z_center = getCenter(n_c[i].elt[2].slice(x));
            Number n_c_z_center_2 = n_c_z_center * n_c_z_center;
            // get radii of independent generators and their squares
            TYPE n_c_x_radius = getRadius(n_c[i].elt[0].slice(x).independent)
            TYPE n_c_x_radius_2 = n_c_x_radius * n_c_x_radius;
            TYPE n_c_y_radius = getRadius(n_c[i].elt[1].slice(x).independent)
            TYPE n_c_y_radius_2 = n_c_y_radius * n_c_y_radius;
            TYPE n_c_z_radius = getRadius(n_c[i].elt[2].slice(x).independent)
            TYPE n_c_z_radius_2 = n_c_z_radius * n_c_z_radius;

            //     separation constraint: -inf < -1*f_c_z < 0
            // storing separation constraint value, not sure what index to use here?
            // question: not sure what to do with the radius here?
            // g[i] = -1*getCenter(f_c_z);
            g[i] = -1*f_c_z_center;

            //     slipping constraint: -inf < f_c_x*f_c_x + f_c_y*f_c_y - u_s^2*f_c_z*f_c_z < 0
            // need to write this constraint differently than matlab implementation as we need to 
            // slice before doing multiplications of PZs in order to avoid memory problems.
            // getCenter, getRadius (for v_norm in armour_main.cpp : getRadius(u_nom[t_ind * NUM_FACTORS + i].independent))

            // if ( getCenter(f_c[i].elt[0].slice(x)) >= 0) && (getCenter(f_c[i].elt[1].slice(x)) >= 0) && (getCenter(f_c[i].elt[2].slice(x)) >= 0 ){
            //     g[i+NUM_TIME_STEPS] = (getCenter(f_c[i].elt[0].slice(x)))_2 + 2*getRadius(f_c[i].elt[0].slice(x).independent)*getCenter(f_c[i].elt[0].slice(x)) + (getRadius(f_c[i].elt[0].slice(x).independent))_2 + (getCenter(f_c[i].elt[1].slice(x)))_2 + 2*getRadius(f_c[i].elt[1].slice(x).independent)*getCenter(f_c[i].elt[1].slice(x)) + (getRadius(f_c[i].elt[1].slice(x).independent))_2 - u_s_2*((getCenter(f_c[i].elt[2].slice(x)))_2 - 2*getRadius(f_c[i].elt[2].slice(x).independent)*getCenter(f_c[i].elt[2].slice(x)) - (getRadius(f_c[i].elt[2].slice(x).independent))_2);
            
            // check the signs of the centers of the force zonotopes
            // condition 1: all positive
            if ( (f_c_x_center >= 0) && (f_c_y_center >= 0) && (f_c_z_center >= 0) ){
                // Note: double check that the center/radius is a number that can be squared
                g[i+NUM_TIME_STEPS] = f_c_x_center_2 + 2*f_c_x_radius*f_c_x_center + f_c_x_radius_2 + f_c_y_center_2 + 2*f_c_y_radius*f_c_y_center + f_c_y_radius_2 - u_s*u_s * ( f_c_z_center_2 - 2*f_c_z_radius*f_c_z_center - f_c_z_radius_2); // checked signs
            }
            // condition 2: y negative
            else if ( (f_c_x_center >= 0) && (f_c_y_center <= 0) && (f_c_z_center >= 0) ) {
                g[i+NUM_TIME_STEPS] = f_c_x_center_2 + 2*f_c_x_radius*f_c_x_center + f_c_x_radius_2 + f_c_y_center_2 - 2*f_c_y_radius*f_c_y_center + f_c_y_radius_2 - u_s*u_s * ( f_c_z_center_2 - 2*f_c_z_radius*f_c_z_center - f_c_z_radius_2); // checked signs
            }
            // condition 3: z negative
            else if ( (f_c_x_center >= 0) && (f_c_y_center >= 0) && (f_c_z_center <= 0) ) {
                g[i+NUM_TIME_STEPS] = f_c_x_center_2 + 2*f_c_x_radius*f_c_x_center + f_c_x_radius_2 + f_c_y_center_2 + 2*f_c_y_radius*f_c_y_center + f_c_y_radius_2 - u_s*u_s * ( f_c_z_center_2 + 2*f_c_z_radius*f_c_z_center - f_c_z_radius_2); // checked signs
            }
            // condition 4: y and z negative
            else if ( (f_c_x_center >= 0) && (f_c_y_center <= 0) && (f_c_z_center <= 0) ) {
                g[i+NUM_TIME_STEPS] = f_c_x_center_2 + 2*f_c_x_radius*f_c_x_center + f_c_x_radius_2 + f_c_y_center_2 - 2*f_c_y_radius*f_c_y_center + f_c_y_radius_2 - u_s*u_s * ( f_c_z_center_2 + 2*f_c_z_radius*f_c_z_center - f_c_z_radius_2); // checked signs
            }
            // condition 5: x negative
            else if ( (f_c_x_center <= 0) && (f_c_y_center >= 0) && (f_c_z_center >= 0) ) {
                g[i+NUM_TIME_STEPS] = f_c_x_center_2 - 2*f_c_x_radius*f_c_x_center + f_c_x_radius_2 + f_c_y_center_2 + 2*f_c_y_radius*f_c_y_center + f_c_y_radius_2 - u_s*u_s * ( f_c_z_center_2 - 2*f_c_z_radius*f_c_z_center - f_c_z_radius_2); // checked signs
            }
            // condition 6: x and y negative
            else if ( (f_c_x_center <= 0) && (f_c_y_center <= 0) && (f_c_z_center >= 0) ) {
                g[i+NUM_TIME_STEPS] = f_c_x_center_2 - 2*f_c_x_radius*f_c_x_center + f_c_x_radius_2 + f_c_y_center_2 - 2*f_c_y_radius*f_c_y_center + f_c_y_radius_2 - u_s*u_s * ( f_c_z_center_2 - 2*f_c_z_radius*f_c_z_center - f_c_z_radius_2); // checked signs
            }
            // condition 7: x and z negative
            else if ( (f_c_x_center <= 0) && (f_c_y_center >= 0) && (f_c_z_center <= 0) ) {
                g[i+NUM_TIME_STEPS] = f_c_x_center_2 - 2*f_c_x_radius*f_c_x_center + f_c_x_radius_2 + f_c_y_center_2 + 2*f_c_y_radius*f_c_y_center + f_c_y_radius_2 - u_s*u_s * ( f_c_z_center_2 + 2*f_c_z_radius*f_c_z_center - f_c_z_radius_2); // checked signs
            }
            // condition 8: x and y and z negative
            else if ( (f_c_x_center <= 0) && (f_c_y_center <= 0) && (f_c_z_center <= 0) ) {
                g[i+NUM_TIME_STEPS] = f_c_x_center_2 - 2*f_c_x_radius*f_c_x_center + f_c_x_radius_2 + f_c_y_center_2 - 2*f_c_y_radius*f_c_y_center + f_c_y_radius_2 - u_s*u_s * ( f_c_z_center_2 + 2*f_c_z_radius*f_c_z_center - f_c_z_radius_2); // checked signs
            }

            //     tipping constraint: ZMP_top_x*ZMP_top_x + ZMP_top_y*ZMP_top_y - surf_rad*ZMP_bottom*ZMP_bottom < 0

            // compute the numerator of the ZMP point equation
            vecPZsparse ZMP_top = cross([0;0;1],n_c[i]);
            // extract the x, y and z components, slice by the parameters, then get the centers and radii of independent generators and their squares
            // x-component
            PZsparse ZMP_top_x = ZMP_top.elt[0].slice(x);
            ZMP_top_x_center = getCenter(ZMP_top_x);
            ZMP_top_x_center_2 = ZMP_top_x_center*ZMP_top_x_center;
            ZMP_top_x_radius = getRadius(ZMP_top_x.independent);
            ZMP_top_x_radius_2 = ZMP_top_x_radius*ZMP_top_x_radius;
            // y-component
            ZMP_top_y = ZMP_top.elt[1].slice(x);
            ZMP_top_y_center = getCenter(ZMP_top_y);
            ZMP_top_y_center_2 = ZMP_top_y_center*ZMP_top_y_center;
            ZMP_top_y_radius = getRadius(ZMP_top_y.independent);
            ZMP_top_y_radius_2 = ZMP_top_y_radius*ZMP_top_y_radius;
            // z-component
            ZMP_top_z = ZMP_top.elt[2].slice(x); // use for debugging, check this is always equal to zero
            
            // compute the denominator of the ZMP point equation
            // note that if the normal vector corresponds to the body frame z-axis, n=[0;0;1] and the dot product of that normal vector
            // with the moment results in the z-component of the moment.
            // question: is the moment about the contact point like in Matlab RNEA? need to verify this
            ZMP_bottom = n_c[i].elt[2].slice(x);
            // extract center and radius of independent generators and their squares
            ZMP_bottom_center = getCenter(ZMP_bottom);
            ZMP_bottom_center_2 = ZMP_bottom_center*ZMP_bottom_center;
            ZMP_bottom_radius = getRadius(ZMP_bottom.independent);
            ZMP_bottom_radius_2 = ZMP_bottom_radius*ZMP_bottom_radius;
            
            // check the signs of the centers of the force zonotopes in order to form the constraint
            // condition 1: all positive
            if ( (ZMP_top_x_center >= 0) && (ZMP_top_y_center >= 0) && (ZMP_bottom_center >= 0) ){
                // Note: double check that the center/radius is a number that can be squared
                g[i+2*NUM_TIME_STEPS] = ZMP_top_x_center_2 + 2*ZMP_top_x_radius*ZMP_top_x_center + ZMP_top_x_radius_2 + ZMP_top_y_center_2 + 2*ZMP_top_y_radius*ZMP_top_y_center + ZMP_top_y_radius_2 - surf_rad_2 * ( ZMP_bottom_center_2 - 2*ZMP_bottom_radius*ZMP_bottom_center - ZMP_bottom_radius_2); // checked signs
            }
            // condition 2: y negative
            else if ( (ZMP_top_x_center >= 0) && (ZMP_top_y_center <= 0) && (ZMP_bottom_center >= 0) ) {
                g[i+2*NUM_TIME_STEPS] = ZMP_top_x_center_2 + 2*ZMP_top_x_radius*ZMP_top_x_center + ZMP_top_x_radius_2 + ZMP_top_y_center_2 - 2*ZMP_top_y_radius*ZMP_top_y_center + ZMP_top_y_radius_2 - surf_rad_2 * ( ZMP_bottom_center_2 - 2*ZMP_bottom_radius*ZMP_bottom_center - ZMP_bottom_radius_2); // checked signs
            }
            // condition 3: z negative
            else if ( (ZMP_top_x_center >= 0) && (ZMP_top_y_center >= 0) && (ZMP_bottom_center <= 0) ) {
                g[i+2*NUM_TIME_STEPS] = ZMP_top_x_center_2 + 2*ZMP_top_x_radius*ZMP_top_x_center + ZMP_top_x_radius_2 + ZMP_top_y_center_2 + 2*ZMP_top_y_radius*ZMP_top_y_center + ZMP_top_y_radius_2 - surf_rad_2 * ( ZMP_bottom_center_2 + 2*ZMP_bottom_radius*ZMP_bottom_center - ZMP_bottom_radius_2); // checked signs
            }
            // condition 4: y and z negative
            else if ( (ZMP_top_x_center >= 0) && (ZMP_top_y_center <= 0) && (ZMP_bottom_center <= 0) ) {
                g[i+2*NUM_TIME_STEPS] = ZMP_top_x_center_2 + 2*ZMP_top_x_radius*ZMP_top_x_center + ZMP_top_x_radius_2 + ZMP_top_y_center_2 - 2*ZMP_top_y_radius*ZMP_top_y_center + ZMP_top_y_radius_2 - surf_rad_2 * ( ZMP_bottom_center_2 + 2*ZMP_bottom_radius*ZMP_bottom_center - ZMP_bottom_radius_2); // checked signs
            }
            // condition 5: x negative
            else if ( (ZMP_top_x_center <= 0) && (ZMP_top_y_center >= 0) && (ZMP_bottom_center >= 0) ) {
                g[i+2*NUM_TIME_STEPS] = ZMP_top_x_center_2 - 2*ZMP_top_x_radius*ZMP_top_x_center + ZMP_top_x_radius_2 + ZMP_top_y_center_2 + 2*ZMP_top_y_radius*ZMP_top_y_center + ZMP_top_y_radius_2 - surf_rad_2 * ( ZMP_bottom_center_2 - 2*ZMP_bottom_radius*ZMP_bottom_center - ZMP_bottom_radius_2); // checked signs
            }
            // condition 6: x and y negative
            else if ( (ZMP_top_x_center <= 0) && (ZMP_top_y_center <= 0) && (ZMP_bottom_center >= 0) ) {
                g[i+2*NUM_TIME_STEPS] = ZMP_top_x_center_2 - 2*ZMP_top_x_radius*ZMP_top_x_center + ZMP_top_x_radius_2 + ZMP_top_y_center_2 - 2*ZMP_top_y_radius*ZMP_top_y_center + ZMP_top_y_radius_2 - surf_rad_2 * ( ZMP_bottom_center_2 - 2*ZMP_bottom_radius*ZMP_bottom_center - ZMP_bottom_radius_2); // checked signs
            }
            // condition 7: x and z negative
            else if ( (ZMP_top_x_center <= 0) && (ZMP_top_y_center >= 0) && (ZMP_bottom_center <= 0) ) {
                g[i+2*NUM_TIME_STEPS] = ZMP_top_x_center_2 - 2*ZMP_top_x_radius*ZMP_top_x_center + ZMP_top_x_radius_2 + ZMP_top_y_center_2 + 2*ZMP_top_y_radius*ZMP_top_y_center + ZMP_top_y_radius_2 - surf_rad_2 * ( ZMP_bottom_center_2 + 2*ZMP_bottom_radius*ZMP_bottom_center - ZMP_bottom_radius_2); // checked signs
            }
            // condition 8: x and y and z negative
            else if ( (ZMP_top_x_center <= 0) && (ZMP_top_y_center <= 0) && (ZMP_bottom_center <= 0) ) {
                g[i+2*NUM_TIME_STEPS] = ZMP_top_x_center_2 - 2*ZMP_top_x_radius*ZMP_top_x_center + ZMP_top_x_radius_2 + ZMP_top_y_center_2 - 2*ZMP_top_y_radius*ZMP_top_y_center + ZMP_top_y_radius_2 - surf_rad_2 * ( ZMP_bottom_center_2 + 2*ZMP_bottom_radius*ZMP_bottom_center - ZMP_bottom_radius_2); // checked signs
            }

        }
        
    }

    // Part 4. check collision between joint position reachable set and obstacles (in gpu)
    obstacles->linkFRSConstraints(checkJointsPosition, nullptr, g + (NUM_TIME_STEPS * NUM_FACTORS + NUM_TIME_STEPS*3), nullptr);

    // Part 5. (position & velocity) state limit constraints
    desired_trajectory->returnJointPositionExtremum(g + (NUM_TIME_STEPS * NUM_FACTORS + NUM_TIME_STEPS*3) + (NUM_FACTORS - 1) * NUM_TIME_STEPS * obstacles->num_obstacles, x);
    desired_trajectory->returnJointVelocityExtremum(g + (NUM_TIME_STEPS * NUM_FACTORS + NUM_TIME_STEPS*3) + (NUM_FACTORS - 1) * NUM_TIME_STEPS * obstacles->num_obstacles + NUM_FACTORS * 2, x);

    return true;
}
// [TNLP_eval_g]


// [TNLP_eval_jac_g]
// return the structure or values of the Jacobian
bool armtd_NLP::eval_jac_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Index         nele_jac,
   Index*        iRow,
   Index*        jCol,
   Number*       values
)
{
    if(n != NUM_FACTORS){
        WARNING_PRINT("*** Error wrong value of n in eval_g!");
    }
    if(m != constraint_number){
        WARNING_PRINT("*** Error wrong value of m in eval_g!");
    }
        
    if( values == NULL ) {
       // return the structure of the Jacobian
       // this particular Jacobian is dense
        for(Index i = 0; i < m; i++){
            for(Index j = 0; j < n; j++){
                iRow[i * n + j] = i;
                jCol[i * n + j] = j;
            }
        }
    }
    else {
        Index i;
        #pragma omp parallel for private(i) schedule(static, (NUM_TIME_STEPS * NUM_FACTORS + NUM_TIME_STEPS)/ NUM_THREADS)
        for(i = 0; i < (NUM_TIME_STEPS * NUM_FACTORS + NUM_TIME_STEPS); i++) { // not *3 for the contact constraints since they're being done at the same time.

            if (i < NUM_TIME_STEPS * NUM_FACTORS) {
            // Part 1. slice the control input PZ to get the center of the input bound, 
            //         while the radius of the input bound has already been incorporated in ipopt constraint bounds
            control_input[i].slice(values + i * NUM_FACTORS, x);

            // Part 2. slice the forward kinematics PZ to get the center of the joint position bound, 
            //         while the radius of the joint position bound has already been considered in Obstacles class
            checkJointsPosition[i * 3    ] = getCenter(joint_position[i * 3    ].slice(x));
            checkJointsPosition[i * 3 + 1] = getCenter(joint_position[i * 3 + 1].slice(x));
            checkJointsPosition[i * 3 + 2] = getCenter(joint_position[i * 3 + 2].slice(x));
            joint_position[i * 3    ].slice(dk_checkJointsPosition + (i * 3    ) * NUM_FACTORS, x);
            joint_position[i * 3 + 1].slice(dk_checkJointsPosition + (i * 3 + 1) * NUM_FACTORS, x);
            joint_position[i * 3 + 2].slice(dk_checkJointsPosition + (i * 3 + 2) * NUM_FACTORS, x);
            }
            else {
                // Part 3. Force constraints
                // how to properly index this part?

                // .slice(##, x) fills the pointer ## with the gradient, so for one instance it fills up by NUM_FACTORS


                // separation constraint
                // .slice(values + i*NUM_FACTORS, x)
                // should be able to use .slice(##, x) here
                f_c[i].elt[2].slice(values + i*NUM_FACTORS, x)

                // slipping constraint
                // .slice(values + (NUM_FACTORS*NUM_TIME_STEPS + NUM_TIME_STEPS)*NUM_FACTORS, x)
                values[(i+NUM_TIME_STEPS)*NUM_FACTORS] = ;
                // create function for return partial of a passed in index? function(g[i])
                for (int j=0; j<NUM_FACTORS; j++) {
                    for (int k=0; k<NUM_FACTORS; k++) {
                        if (j==k) {
                            if (degreeArray[j]==0) {

                            }
                            else {
                                
                            }
                        }
                        else {

                        }
                    }
                }

                //tipping constraint
                // .slice(, x)
                values[(i+2*NUM_TIME_STEPS)*NUM_FACTORS] = ; // need to fill in all the factors here
                
                
            }
            

        }

        // Part 4. check collision between joint position reachable set and obstacles (in gpu)
        obstacles->linkFRSConstraints(checkJointsPosition, dk_checkJointsPosition, nullptr, values + NUM_TIME_STEPS * NUM_FACTORS * NUM_FACTORS); // needs to be updated

        // Part 5. (position & velocity) state limit constraints
        desired_trajectory->returnJointPositionExtremumGradient(values + (NUM_TIME_STEPS * NUM_FACTORS + (NUM_FACTORS - 1) * NUM_TIME_STEPS * obstacles->num_obstacles) * NUM_FACTORS, x); // needs to be updated
        desired_trajectory->returnJointVelocityExtremumGradient(values + (NUM_TIME_STEPS * NUM_FACTORS + (NUM_FACTORS - 1) * NUM_TIME_STEPS * obstacles->num_obstacles + NUM_FACTORS * 2) * NUM_FACTORS, x); // needs to be updated
    }

    return true;
}
// [TNLP_eval_jac_g]


// void Number::contactConstraintPartial(TYPE* gradient, const Number contactEquation) {
//     // inputs: 
//     //    pointer to storage location?
//     //    constraint to take derivative of? (slip, tip? then have if statement for each)
//     //    pointer to f_c and n_c?

//     if ( contactEquation == 1 ) {
//         // take gradient of slip constraint equation wrt trajectory parameters
        
//     }
//     else if ( contactEquation == 2 ) {

//     }


// }


// [TNLP_eval_h]
//return the structure or values of the Hessian
bool armtd_NLP::eval_h(
   Index         n,
   const Number* x,
   bool          new_x,
   Number        obj_factor,
   Index         m,
   const Number* lambda,
   bool          new_lambda,
   Index         nele_hess,
   Index*        iRow,
   Index*        jCol,
   Number*       values
)
{
    return false;
}
// [TNLP_eval_h]


// [TNLP_finalize_solution]
void armtd_NLP::finalize_solution(
    SolverReturn               status,
    Index                      n,
    const Number*              x,
    const Number*              z_L,
    const Number*              z_U,
    Index                      m,
    const Number*              g,
    const Number*              lambda,
    Number                     obj_value,
    const IpoptData*           ip_data,
    IpoptCalculatedQuantities* ip_cq
)
{
    // here is where we would store the solution to variables, or write to a file, etc
    // so we could use the solution.

    // store the solution
    for( Index i = 0; i < n; i++ ) {
        solution[i] = (TYPE)x[i];
    }

    // check constraint violation manually for Maximum_CpuTime_Exceeded case
    memcpy(g_copy, g, m * sizeof(Number));

    feasible = true;

    // control input constraints
    for( Index i = 0; i < NUM_TIME_STEPS; i++ ) {
        for( Index j = 0; j < NUM_FACTORS; j++ ) {
            if (g_copy[i * NUM_FACTORS + j] < -torque_limits[j] + v_norm[i * NUM_FACTORS + j] - TORQUE_INPUT_CONSTRAINT_VIOLATION_THRESHOLD || 
                g_copy[i * NUM_FACTORS + j] > torque_limits[j] - v_norm[i * NUM_FACTORS + j] + TORQUE_INPUT_CONSTRAINT_VIOLATION_THRESHOLD) {
                feasible = false;
                cout << "        CUDA & C++: Ipopt: Control torque of joint " << j << " at time interval " << i << " exceeds limit!\n";
                cout << "                        value: " << g_copy[i * NUM_FACTORS + j] << "\n";
                cout << "                        range: [ " << -torque_limits[j] + v_norm[i * NUM_FACTORS + j] << ", "
                                                            << torque_limits[j] - v_norm[i * NUM_FACTORS + j] << " ]\n";
                return;
            }
        }
    }    

    // NOTE: need to add force constraints here and adjust the indices offsets after this
    Index offset = NUM_FACTORS * NUM_TIME_STEPS;
    for( Index i = offset; i < NUM_TIME_STEPS; i++) {
        // separation constraint
        if(g_copy[i] > SEPARATION_CONSTRAINT_VIOLATION_THRESHOLD) {
            feasible = false;
            double t_violation = i - offset;
            cout << "        CUDA & C++: Ipopt: Separation constraint violated at time interval: " << t_violation << " \n";
        }
        // slipping constraint
        if(g_copy[i+NUM_TIME_STEPS] > SLIPPING_CONSTRAINT_VIOLATION_THRESHOLD){
            feasible = false;
            double t_violation = i - offset;
            cout << "        CUDA & C++: Ipopt: Slipping constraint violated at time interval: " << t_violation << " \n";
        }
        // tipping constraint
        if(g_copy[i+2*NUM_TIME_STEPS] > TIPPING_CONSTRAINT_VIOLATION_THRESHOLD){
            feasible = false;
            double t_violation = i - offset;
            cout << "        CUDA & C++: Ipopt: Tipping constraint violated at time interval: " << t_violation << " \n";
        }
    }
    offset +=  NUM_TIME_STEPS*3;

    // collision avoidance constraints
    // Index offset = NUM_FACTORS * NUM_TIME_STEPS + NUM_TIME_STEPS*3;
    for( Index i = 0; i < NUM_FACTORS - 1; i++ ) {
        for( Index j = 0; j < NUM_TIME_STEPS; j++ ) {
            for( Index h = 0; h < obstacles->num_obstacles; h++ ) {
                if (g_copy[(i * NUM_TIME_STEPS + j) * obstacles->num_obstacles + h + offset] > COLLISION_AVOIDANCE_CONSTRAINT_VIOLATION_THRESHOLD) {
                    feasible = false;
                    cout << "        CUDA & C++: Ipopt: Collision between link " << i + 1 << " and obstacle " << h << " at time interval " << j << "!\n";
                    cout << "                        value: " << g_copy[(i * NUM_TIME_STEPS + j) * obstacles->num_obstacles + h + offset] << "\n";
                    return;
                }
            }
        }
    }
    offset += (NUM_FACTORS - 1) * NUM_TIME_STEPS * obstacles->num_obstacles;

    

    // state limit constraints
    //     minimum joint position
    for( Index i = offset; i < offset + NUM_FACTORS; i++ ) {
        if (g_copy[i] < state_limits_lb[i - offset] + qe || g_copy[i] > state_limits_ub[i - offset] - qe) {
            feasible = false;
            cout << "        CUDA & C++: Ipopt: joint " << i - offset << " exceeds position limit when it reaches minimum!\n";
            cout << "                        value: " << g_copy[i] << "\n";
            cout << "                        range: [ " << state_limits_lb[i - offset] + qe << ", "
                                                        << state_limits_ub[i - offset] - qe << " ]\n";
            return;
        }
    }
    offset += NUM_FACTORS;

    //     maximum joint position
    for( Index i = offset; i < offset + NUM_FACTORS; i++ ) {
        if (g_copy[i] < state_limits_lb[i - offset] + qe || g_copy[i] > state_limits_ub[i - offset] - qe) {
            feasible = false;
            cout << "        CUDA & C++: Ipopt: joint " << i - offset << " exceeds position limit when it reaches maximum!\n";
            cout << "                        value: " << g_copy[i] << "\n";
            cout << "                        range: [ " << state_limits_lb[i - offset] + qe << ", "
                                                        << state_limits_ub[i - offset] - qe << " ]\n";
            return;
        }
    }
    offset += NUM_FACTORS;

    //     minimum joint velocity
    for( Index i = offset; i < offset + NUM_FACTORS; i++ ) {
        if (g_copy[i] < -speed_limits[i - offset] + qde || g_copy[i] > speed_limits[i - offset] - qde) {
            feasible = false;
            cout << "        CUDA & C++: Ipopt: joint " << i - offset << " exceeds velocity limit when it reaches minimum!\n";
            cout << "                        value: " << g_copy[i] << "\n";
            cout << "                        range: [ " << -speed_limits[i - offset] + qde << ", "
                                                        << speed_limits[i - offset] - qde << " ]\n";
            return;
        }
    }
    offset += NUM_FACTORS;

    //     maximum joint velocity
    for( Index i = offset; i < offset + NUM_FACTORS; i++ ) {
        if (g_copy[i] < -speed_limits[i - offset] + qde || g_copy[i] > speed_limits[i - offset] - qde) {
            feasible = false;
            cout << "        CUDA & C++: Ipopt: joint " << i - offset << " exceeds velocity limit when it reaches maximum!\n";
            cout << "                        value: " << g_copy[i] << "\n";
            cout << "                        range: [ " << -speed_limits[i - offset] + qde << ", "
                                                        << speed_limits[i - offset] - qde << " ]\n";
            return;
        }
    }
}
// [TNLP_finalize_solution]


#endif
