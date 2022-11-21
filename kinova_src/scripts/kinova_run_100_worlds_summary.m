%% description
% This script gets the statics of the summaries generated from
% kinova_run_100_worlds
%
% Authors: Bohao Zhang (adapted from Shreyas Kousik code)
% Created 25 October 2022

clear; clc;

use_robust_input = true;

save_file_header = 'trial_' ;
file_location = '../simulator_files/testing/saved_worlds/20221023_armtd_withoutgripper' ;
if use_robust_input
    file_location = [file_location, '_robust'];
else
    file_location = [file_location, '_nominal_only'];
end
addpath(file_location);

summary_files = dir([file_location, '/trial_*']);

collision_check = [];
input_check = [];
ultimate_bound_check = [];
joint_limit_check = [];
goal_check = [];
infeasible_check = [];

for i = 1:length(summary_files)
    data = load(summary_files(i).name);
    if data.summary.collision_check
        collision_check = [collision_check, i];
        continue;
    end
    if data.summary.input_check
        input_check = [input_check, i];
        continue;
    end
    if data.summary.ultimate_bound_check
        ultimate_bound_check = [ultimate_bound_check, i];
        continue;
    end
    if data.summary.joint_limit_check
        joint_limit_check = [joint_limit_check, i];
        continue;
    end
    if data.summary.goal_check
        goal_check = [goal_check, i];
        continue;
    end
    infeasible_check = [infeasible_check, i];
end

collision_check
input_check
ultimate_bound_check
joint_limit_check
goal_check
infeasible_check