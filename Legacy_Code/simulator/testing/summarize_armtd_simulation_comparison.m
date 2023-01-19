%% description
% This script iterates through the saved armtd trials and summarizes
% information about each (ie crash check, goal check, planning time, etc.)
%
% Author: Patrick Holmes
% Created 06 December 2019

clear; clc;

trial_file_header = 'trial_scene_' ;
% trial_file_folder = '../armtd-dev/simulator_files/testing/trial_data/20200127_orig_armtd/' ;
% trial_file_folder = './trial_data/20211217/' ;
% trial_file_folder = './trial_data/20211122_new_cuda_UB_dynamics_spin/' ; 
% trial_file_folder = '../armtd-dev/simulator_files/testing/trial_data/20211122_new_cuda_UB_dynamics_spin/' ;
% trial_file_folder = '../armtd-dev/simulator_files/testing/trial_data/20211122_new_cuda_UB_dynamics_quarks/' ;
% trial_file_folder = '../armtd-dev/simulator_files/testing/trial_data/20211122_new_cuda_UB/' ;
% trial_file_folder = '../armtd-dev/simulator_files/testing/trial_data/20211122_new_robust_input_no_obstacles_20220127/';
% trial_file_folder = '../armtd-dev/simulator_files/testing/trial_data/20220602_fetch_dumbbell_input_constraints_nominal_only/';
% trial_file_folder = '../armtd-dev/simulator_files/testing/trial_data/20220602_fetch_dumbbell_input_constraints_longer_timeout_robust/';
% trial_file_folder = '../armtd-dev/simulator_files/testing/trial_data/20220602_fetch_dumbbell/';
% trial_file_folder = '../armtd-dev/simulator_files/testing/trial_data/20220602_fetch_dumbbell_no_input_constraints_robust/';
% trial_file_folder = '../armtd-dev/simulator_files/testing/trial_data/20220602_fetch_dumbbell_no_input_constraints_true_params_for_robust/';
trial_file_location = sprintf('%s*%s*', trial_file_folder, trial_file_header);
trial_file_list = dir(trial_file_location);

results = struct();
results.scene = {};
results.crash = [];
results.inputs = [];
results.ultimate_bound = [];
results.joint_limit = [];
results.goal = [];
results.avg_planning_time = [];
results.time_to_goal = [];
results.timed_out = [];

for idx = 1:length(trial_file_list)
    trial_filename = trial_file_list(idx).name;
    mytrial = load([trial_file_folder trial_filename]);
    summary = mytrial.summary;
    
    results.scene{end+1, 1} = mytrial.world_filename;
    results.crash(end+1, 1) = summary.collision_check;
    results.inputs(end+1, 1) = summary.input_check;
    results.ultimate_bound(end+1, 1) = summary.ultimate_bound_check;
    results.joint_limit(end+1, 1) = summary.joint_limit_check;
    results.goal(end+1, 1) = summary.goal_check;
    results.timed_out(end+1, 1) = summary.total_real_time > 172800;
    if summary.goal_check
       results.avg_planning_time(end+1, 1) = mean(summary.planning_time(~isnan(summary.planning_time)));
    else
       results.avg_planning_time(end+1, 1) = nan; 
    end
    if summary.goal_check
       results.time_to_goal(end+1, 1) = summary.total_simulated_time(end);
    else
       results.time_to_goal(end+1, 1) = nan; 
    end
end

disp(struct2table(results));

disp('Crashed?');
disp([num2str(sum(results.crash)) ' / ' num2str(length(results.crash))]);
disp('Exceeded input bounds?');
disp([num2str(sum(results.inputs)) ' / ' num2str(length(results.inputs))]);
disp('Exceeded ultimate bound?');
disp([num2str(sum(results.ultimate_bound)) ' / ' num2str(length(results.ultimate_bound))]);
disp('Exceeded joint limits?');
disp([num2str(sum(results.joint_limit)) ' / ' num2str(length(results.joint_limit))]);
disp('Reached Goal?');
disp([num2str(sum(results.goal)) ' / ' num2str(length(results.goal))]);

disp('Timed out?');
% disp([num2str(sum(results.timed_out)) ' / ' num2str(length(results.timed_out))]);
% disp('Avg. planning time (when goal reached)');
% disp(num2str(nanmean(results.avg_planning_time)));
% disp('Avg. time to goal (when goal reached):');
% disp(num2str(nanmean(results.time_to_goal)));
