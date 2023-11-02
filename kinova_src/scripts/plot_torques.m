%%
function plot_torques(input, input_centers, input_radii, time, trajectory)

    iterations = size(input_centers, 1);
    JOINTS = size(input, 1);
    disp(JOINTS)

    disp(iterations)
    
    % Reshape the input to have 100 timesteps for each iteration
    trajectories = reshape(input(:, 1:end-1), [JOINTS, 100, iterations]);
    time = reshape(time(1:end-1), [100, iterations]);
    if trajectory == -1
    else
        trajectories = trajectories(:, :, trajectory);
        time = time(:, trajectory);
        iterations = 1;
    end
    trajectories = permute(trajectories, [3, 1, 2]);  % (iterations, JOINTS, 100)
    
    % Prepare the new time scale for the centers/radii
    old_time_scale = linspace(0, 1, 128);
    new_time_scale = linspace(0, 1, 100);
    
    % Scale the torque_centers and torque_radii
    scaled_torque_centers = zeros(iterations, JOINTS, 100);
    scaled_torque_radii = zeros(iterations, JOINTS, 100);
    
    for j = 1:JOINTS
        for i = 1:iterations
            scaled_torque_centers(i, j, :) = interp1(old_time_scale, squeeze(input_centers(i, j, :))', new_time_scale, 'linear', 'extrap');
            scaled_torque_radii(i, j, :) = interp1(old_time_scale, squeeze(input_radii(i, j, :))', new_time_scale, 'linear', 'extrap');
        end
    end
    
    % Loop through each iteration to plot
    for i = 1:iterations
        figure;  % Creates a new figure for each iteration
        hold on; % Hold on to plot multiple datasets
        x = squeeze(time(:, i)');
        if trajectory == -1
            sgtitle("Plots for trajectory: " + i)
        else
            sgtitle("Plots for trajectory: " + trajectory)
        end
        
        
        % Loop to calculate and plot the bounds and input for each joint
        for j = 1:JOINTS
            % Select the j-th subplot
            subplot(ceil(JOINTS/2), 2, j);
            
            % Calculate bounds for this joint
            upper_bound = squeeze(scaled_torque_centers(i, j, :) + scaled_torque_radii(i, j, :))';
            lower_bound = squeeze(scaled_torque_centers(i, j, :) - scaled_torque_radii(i, j, :))';
            x2 = [x, fliplr(x)];
            inBetween = [lower_bound, fliplr(upper_bound)];
            
            % Shade the area between bounds
            fill(x2, inBetween, 'k', 'facealpha', 0.1); % Shading
            hold on; % Hold on to plot the input on the same subplot
            
            % Plot the input trajectory on top
            plot(x, squeeze(trajectories(i, j, :)), 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]); % Input trajectory
            
            % Labels and titles for each subplot
            xlabel('Time');
            title(['Joint ' num2str(j)]);
            
            hold off; % Release hold for the next subplot
        end
        
    end
end