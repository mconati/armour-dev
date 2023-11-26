
function plot_bounded_trajectories(plt, A, T, Y, bounds, traj_len, braking_frac, joint)
    axes(plt); 
    hold on;    
    iterations = (length(A.time)-1)*2/traj_len;
    input_shape = size(Y);
    joints = input_shape(1);
    if joint<joints
        Y = Y(joint, :);
        joints = 1;
    end
    [moving, braking] = parse_full_trajectories(T, Y, iterations, traj_len, braking_frac);
    moving_shape = size(moving);
    moving_time = A.full_time(1:moving_shape(2));
    %disp(size(braking))
    plot(plt, moving_time, moving(:, :), color='b')
    hold on;
    offset = 1;
    for i = 1:1:iterations
        start_index = (i)*braking_frac*traj_len+1;
        end_index = (i)*braking_frac*traj_len+1+(traj_len-braking_frac*traj_len);
        braking_time = A.full_time(start_index:end_index-1);
        plot(plt, braking_time, braking(offset:offset+joints-1, :), color='r')
        offset = offset+joints;
        hold on;
    end
    for i=1:1:length(bounds)/2
        y_min = bounds(2*i-1); % Lower bound for the i-th subplot
        y_max = bounds(2*i);   % Upper bound for the i-th subplot
        plot(plt, T/2, repmat(y_min, 1, length(T)), 'k');
        plot(plt, T/2, repmat(y_max, 1, length(T)), 'g');
        hold on;
    end
    
function [moving, braking] = parse_full_trajectories(T, Y, iterations, traj_len, braking_frac)
    moving = [];
    braking = [];
    for i = 1:1:iterations
        start_index = (i-1)*traj_len+2;
        end_index = (i)*traj_len+2;
        mid_index = braking_frac * traj_len + start_index;
        moving = [moving, Y(:, start_index:mid_index-1)];
        braking = [braking; Y(:, mid_index:end_index-1)];
    end
end
end




