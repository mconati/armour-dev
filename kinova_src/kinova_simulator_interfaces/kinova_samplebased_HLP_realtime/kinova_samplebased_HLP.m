classdef kinova_samplebased_HLP < robot_arm_graph_planner_HLP
    properties
        sample_nodes = [];
        start_goal_distance_threshold = norm(pi / 48 * ones(1,7));
    end

    methods
        function HLP = kinova_samplebased_HLP(varargin)
            HLP@robot_arm_graph_planner_HLP(varargin{:}) ;
%             HLP.sample_nodes = load('millionNodes.mat');
            HLP.sample_nodes = load('../kinova_simulator_interfaces/kinova_samplebased_HLP_realtime/uniformNodes.mat');
        end

        function HLP = generatePath(HLP, obstacles, start, goal)
            % obstacle number hardcoded as 10
            Zs = [];
            for i = 1:10
                Zs = [Zs; obstacles{i}.Z'];
            end
            
            % write obstacle info to file as the input of the CUDA collision checker
            writematrix(Zs, '../kinova_simulator_interfaces/kinova_samplebased_HLP_realtime/obstacles.csv', 'Delimiter', ' ');
            
            % call collision checker in CUDA
            system('./../kinova_simulator_interfaces/kinova_samplebased_HLP_realtime/collision_checker');
            
            adj_matrix_sparse_data = readmatrix('../kinova_simulator_interfaces/kinova_samplebased_HLP_realtime/collision_free_adj_matrix.csv');
            adj_matrix_sparse = sparse(adj_matrix_sparse_data(:,1)+1, ...
                                       adj_matrix_sparse_data(:,2)+1, ...
                                       adj_matrix_sparse_data(:,3), ...
                                       size(HLP.sample_nodes.q_valid_list,2), size(HLP.sample_nodes.q_valid_list,2));
            G = graph(adj_matrix_sparse, 'lower');
            
            [bins, binsize] = conncomp(G);
            [~, max_id] = max(binsize);
            G_maxconn = subgraph(G, bins == max_id);
            
            q_subgraph = HLP.sample_nodes.q_valid_list(:, bins == max_id);
            difference_to_goal = vecnorm(wrapToPi(goal - q_subgraph));
            difference_to_start = vecnorm(wrapToPi(start - q_subgraph));
            [end_diff, end_idx] = min(difference_to_goal);
            [start_diff, start_idx] = min(difference_to_start);

            if end_diff > HLP.start_goal_distance_threshold
                fprintf('    HLP: Goal node is far away!!! Distance: %f\n', start_diff)
            end
            if start_diff > HLP.start_goal_distance_threshold
                fprintf('    HLP: Start node is far away!!! Distance: %f\n', end_diff)
            end
            
            [path, len] = shortestpath(G_maxconn, start_idx, end_idx);
            
            if isempty(path)
                error('can not find any path!');
            end

            HLP.graph_waypoints = q_subgraph(:,path);
        end

        function plot(HLP)
            % plot the waypoint
            if ~isempty(HLP.current_waypoint_patch_data)
                HLP.plot_data.waypoint_arm_volume = patch(HLP.current_waypoint_patch_data,...
                    'FaceColor','g','FaceAlpha',0.1) ;
            end
        end
    end
end



