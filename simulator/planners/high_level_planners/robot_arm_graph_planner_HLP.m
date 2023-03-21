classdef robot_arm_graph_planner_HLP < high_level_planner
    %% Description
    % High level planner that takes in a list of waypoints that are nodes
    % in a graph. This planner returns a waypoint using a straight line
    % planner to get to the current closest passed in waypoint.
    
    %% properties
    properties

        % arm
        arm_n_states
        arm_n_inputs
        arm_n_links_and_joints
        arm_joint_state_limits
        arm_joint_speed_limits
        arm_joint_input_limits
        arm_joint_state_indices
        arm_joint_speed_indices
        
        % samples
%         nodes
%         adjacency_matrix
%         n_nodes
%         n_nodes_max = 20000 ;
%         all_node_idxs ;
%         sampling_timeout = 0.1 ; % seconds per planning iteration
        
        % path
%         best_path_nodes
%         best_path_node_idxs

        graph
        graph_waypoints
        current_graph_waypoint_index = 0;

    end
    methods
        %%  constructor
        function HLP = robot_arm_graph_planner_HLP(varargin)
            HLP@high_level_planner(varargin{:}) ;
        end
        
        %% setup
        function setup(HLP,agent_info,world_info)
            % get all the necessary arm properties filled in
            HLP.vdisp('Filling in HLP''s arm properties',9)
            HLP = fill_in_arm_properties(HLP,agent_info,false) ;
            
            % get the world goal
            HLP.vdisp('Setting goal',9)
            HLP.goal = world_info.goal ;
        end
        
        %% get waypoint
        function waypoint = get_waypoint(HLP,agent_info,~,lookahead_distance)

            q_cur = agent_info.state(HLP.arm_joint_state_indices, end);

            if isempty(HLP.current_waypoint) % meaning just starting

                % increment current graph waypoint index
%                 graph_waypoint = HLP.graph_waypoints(:,1);
                HLP.current_graph_waypoint_index = 1;

                % call Straight Line Planner function to get waypoint in
                % direction of graph waypoint
                waypoint = get_SLP_to_waypoint(HLP,agent_info,lookahead_distance);

                % update current waypoint (done in the SLP function?)
%                 HLP.current_waypoint = waypoint;
%                 HLP.current_waypoint_index = 1;

            else

                % check how close to current waypoint
                cur_dist = norm(q_cur-HLP.current_waypoint);

                % if too far away, return SLP waypoint based on current
                % graph waypoint
                if cur_dist > lookahead_distance
                    waypoint = get_SLP_to_waypoint(HLP,agent_info,lookahead_distance);

                else % choose next waypoint
                    if HLP.current_graph_waypoint_index == length(HLP.graph_waypoints) % (:,end) % if already at last waypoint, return goal
                        waypoint = HLP.goal;
                        HLP.current_waypoint_index = length(HLP.waypoints);
                    else % return the next waypoint
                        HLP.current_graph_waypoint_index = HLP.current_graph_waypoint_index + 1;
                        waypoint = get_SLP_to_waypoint(HLP,agent_info,lookahead_distance);
%                         HLP.current_waypoint = waypoint;
                    end
                end
            end
        end

        function waypoint = get_SLP_to_waypoint(HLP,agent_info,lookahead_distance)
            if nargin < 3
                lookahead_distance = HLP.default_lookahead_distance;
            end
    
            g = HLP.graph_waypoints(:,HLP.current_graph_waypoint_index);
            z = agent_info.state(HLP.arm_joint_state_indices, end);
            dir_des = g - z;
            dir_des = dir_des./norm(dir_des);
            
            % adjust lookahead distance to not overshoot goal
            dist_to_goal = norm(g - z);
            if lookahead_distance > dist_to_goal
                lookahead_distance = dist_to_goal;
            end
            
            waypoint = lookahead_distance.*dir_des + z;
    
            % update current waypoints
            HLP.current_waypoint = waypoint;
            HLP.waypoints = [HLP.waypoints, waypoint];
            HLP.current_waypoint_index = HLP.current_waypoint_index + 1;
        end


    end
end