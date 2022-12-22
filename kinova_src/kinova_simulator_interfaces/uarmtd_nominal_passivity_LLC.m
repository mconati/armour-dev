classdef uarmtd_nominal_passivity_LLC < robot_arm_LLC
    % utilizes nominal params with RNEA
    % adapted from uarmtd_robust_CBF_LLC, just taking out robust part.
    % but keeping ultimate bound checking
    
    properties

    end
    
    methods
        function LLC = uarmtd_nominal_passivity_LLC(varargin)
            LLC = parse_args(LLC,varargin{:}) ;
        end
        
        %% setup
        function setup(LLC,A)
            % call default setup
            setup@robot_arm_LLC(LLC, A);
        end

        %% info
        function info = get_LLC_info(LLC)
            info = struct();
        end

        %% get control inputs
        function [u, tau, v, true_disturbance, true_V, r] = get_control_inputs(LLC, A, t, z_meas, planner_info)
            % u = LLC.get_control_inputs(A, t_cur,z_cur,planner_info)
            %
            % Given the current time and measured state, and a reference trajectory
            % as a function of time within planner_info, compute the control
            % input that attempts to track the given reference.

            % get actual robot state
            z = z_meas(:);
            q = z(A.joint_state_indices);
            qd = z(A.joint_speed_indices);

            % compute desired trajectory
            [q_des, qd_des, qdd_des] = planner_info.desired_trajectory{end}(t);

            % error terms
            err = q_des - q;
            d_err = qd_des - qd;

            % modified reference trajectory
            qd_ref = qd_des + LLC.Kr * err;
            qdd_ref = qdd_des + LLC.Kr * d_err;

            r = d_err + LLC.Kr*err;

            % nominal controller
            tau = rnea(q, qd, qd_ref, qdd_ref, true, A.params.nominal);

            % input is just nominal input, no robust input.
            u = tau;

            % output more for logging if desired
            if nargout > 3
                    disturbance = rnea(q, qd, qd_ref, qdd_ref, true, A.params.true) - tau;
                    V_tmp = rnea(q, zeros(A.n_states/2, 1), zeros(A.n_states/2, 1), r, false, A.params.true);
                    V = 0.5*r'*V_tmp;
                    true_disturbance = disturbance;
                    true_V = V;
            end
        end
        
        %% helper functions
        function out = ultimate_bound_check(LLC, A, t_start)
            % create time vector for checking
            t_agent = A.time(end);
            t_check = t_start:A.traj_check_time_discretization:t_agent;

            if isempty(t_check) || t_check(end) ~= t_agent
                t_check = [t_check, t_agent] ;
            end

            % get agent state and reference trajectories interpolated to time
            z_agent = match_trajectories(t_check,A.time,A.state) ;
            z_ref_agent = match_trajectories(t_check,A.time,A.reference_state) ;

            % check bound satisfaction
            A.vdisp('Running ultimate bound check!',3);
            out = false;
            for t_idx = 1:length(t_check)
                q = z_agent(A.joint_state_indices, t_idx);
                qd = z_agent(A.joint_speed_indices, t_idx);
                q_ref = z_ref_agent(A.joint_state_indices, t_idx);
                qd_ref = z_ref_agent(A.joint_speed_indices, t_idx);
                for i = 1:length(q)
                    if abs(q(i) - q_ref(i)) > LLC.ultimate_bound_position
                        fprintf('Time %.2f, joint %d position bound exceeded: %.5f vs +-%.5f \n', t_check(t_idx), i, q(i) - q_ref(i), LLC.ultimate_bound_position);
                        out = true;
                    end
                    if abs(qd(i) - qd_ref(i)) > LLC.ultimate_bound_velocity
                        fprintf('Time %.2f, joint %d velocity bound exceeded: %.5f vs +-%.5f \n', t_check(t_idx), i, qd(i) - qd_ref(i), LLC.ultimate_bound_velocity);
                        out = true;
                    end
                end
            end

            if ~out
                LLC.vdisp('No ultimate bound exceeded', 3);
            end
        end
    end
        
end

