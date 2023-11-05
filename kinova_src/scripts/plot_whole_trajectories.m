
function plot_whole_trajectories(A, T, Y, bounds, makefigure)
    buffer_y = [Y(:, 1)];
    buffer_yt = [T(:, 1)];
    buffer_b = [];
    buffer_bt = [];
    t = A.time;
    
    t_stop = A.t_stop; 
    
    % Create a figure and axis
    if makefigure
       figure;
    end
    hold on;
    
    color_idx = 1; % Initialize the color index
    start_time = 0;
    s_T = 1;
    while start_time < max(T)
        
        end_time = start_time + t_stop;
        % Clip the time interval to avoid going beyond the data
        if end_time > max(T)
            end_time = max(T);
        end
        indices = find(T > start_time & T <= end_time);
        % Select the appropriate color
        if color_idx == 1
            color = 'b'; % Blue
            buffer_y = [buffer_y, Y(:, indices)];
            buffer_yt = [buffer_yt, T(:, s_T:(s_T+length(indices))-1)];
            s_T = s_T + length(indices)-1;
        else
            color = 'r'; % Red
            buffer_b = [buffer_b, Y(:, indices)];
            buffer_bt = [buffer_bt, T(:, s_T:(s_T+length(indices))-1)];
            plot(T(:, s_T:(s_T+length(indices)-1)), Y(:, indices), color);
        end
        
        % Find the corresponding indices in the time array
    
        
        % Plot the line segment
    
        
        % Update the color index and start time for the next segment
        color_idx = 3 - color_idx; % Toggle between 1 and 2
        start_time = end_time;
    
    end
        % Plot the bounds 
    if bounds == -1
        title("Trajectory Plotting with Braking")
    else
        for i=1:1:length(bounds)/2
            y_min = bounds(2*i-1); % Lower bound for the i-th subplot
            y_max = bounds(2*i);   % Upper bound for the i-th subplot
            if y_min<-5 || y_max>5
                continue
            end
            plot(T/2, repmat(y_min, 1, length(T)), 'k');
            plot(T/2, repmat(y_max, 1, length(T)), 'g');
        end
        title("Trajectory Plotting with Braking and Bounds")
    end

    plot(buffer_yt, buffer_y, 'b');
    xlabel('Time (seconds)');
    ylabel('Y');
    grid on;
end