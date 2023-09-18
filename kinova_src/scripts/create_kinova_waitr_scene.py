import os
import json
import argparse
from pathlib import Path
import time
import math
from datetime import datetime
from datetime import datetime
import numpy as np
import pandas as pd

import pybullet as p
import pybullet_data

from pyBulletSimRecorder import PyBulletRecorder

def create_static_scene(scene_filename, robot_basename, use_initial=True, use_current=False, use_bounding = True, use_goal=True, use_waypoint=False):
     # get scene information from csv file
    scene_info = load_scene_info(scene_filename)

    # set collision groups
    env_collision_group = 100
    enable_env_collision = 0

    rob_collision_group = 200
    enable_rob_collision = 0
 
    # load env
    env_name, env_id = load_env(scene_filename, env_collision_group, enable_env_collision, use_floor=False, use_mat=False)

    '''TODO: clean this section up'''
    # create robot names
    init_filename = robot_basename + '_initial.urdf'
    curr_filename = robot_basename + '_current.urdf'
    waypoint_filename = robot_basename + '_waypoint.urdf'
    goal_filename = robot_basename + '_goal.urdf'
    bounding_filename = robot_basename + '_bounding.urdf'

    rob_name = []
    rob_id = []

    # load robots
    init_config = list(scene_info.iloc[0])
    goal_config = list(scene_info.iloc[1])

    if use_initial:
        init_id = load_robot(init_filename, init_config, rob_collision_group, enable_rob_collision)
        rob_name.append(init_filename)
        rob_id.append(init_id)
    if use_current:
        curr_id = load_robot(curr_filename, init_config, rob_collision_group, enable_rob_collision)
        rob_name.append(curr_filename)
        rob_id.append(curr_id)
    if use_bounding:
        bounding_id = load_robot(bounding_filename, init_config, rob_collision_group, enable_rob_collision)
        rob_name.append(bounding_filename)
        rob_id.append(bounding_id)
    if use_waypoint:
        waypoint_id = load_robot(waypoint_filename, goal_config, rob_collision_group, enable_rob_collision)
        rob_name.append(waypoint_filename)
        rob_id.append(waypoint_id)
    if use_goal:
        goal_id = load_robot(goal_filename, goal_config, rob_collision_group, enable_rob_collision)
        rob_name.append(goal_filename)
        rob_id.append(goal_id)
    
    # load obstacles
    obs_name, obs_id = load_obstacles(scene_filename, env_collision_group, enable_env_collision)

    # save names and ids
    names = env_name + rob_name + obs_name
    ids = env_id + rob_id + obs_id

    return names, ids

def create_recorder(names, ids):
    '''Save assets for rendering in Blender'''

    assert len(names) == len(ids)

    recorder = PyBulletRecorder()

    for i in range(len(names)):
        recorder.register_object(ids[i], names[i])

    return recorder

def load_scene_info(scene_name):
    return pd.read_csv(scene_name, header=None)

def load_env(scene_filename, collision_group, enable_collision, use_floor=True, use_walls=False, use_mat=True):
    # get scene info TODO: remove this
    scene_info = load_scene_info(scene_filename)

    base_robot_path = 'assets/robots/'
    base_obs_path = 'assets/objects/'

    floor_filename = base_obs_path + 'floor.urdf'
    wall_filename = base_obs_path + 'wall.urdf'
    mat_filename = base_obs_path + 'mat.urdf'

    names, ids = [], []

    if use_floor:
        names.append(floor_filename)
        ids.append(p.loadURDF(floor_filename, useFixedBase=True, basePosition=[0,0,-0.05]))
        p.setCollisionFilterGroupMask(ids[-1],-1, collision_group, enable_collision)
    if use_walls:
        names.append(wall_filename)
        ids.append(p.loadURDF(wall_filename, useFixedBase=True, basePosition=[0,0,-0.05]))
        p.setCollisionFilterGroupMask(ids[-1],-1, collision_group, enable_collision)
    if use_mat:
        names.append(mat_filename)
        ids.append(p.loadURDF(mat_filename, useFixedBase=True, basePosition=[0,0,-0.01]))
        p.setCollisionFilterGroupMask(ids[-1],-1, collision_group, enable_collision)

    return names, ids

def load_obstacles(scene_filename, collision_group, enable_group_collision,offset=None):
    # get scene info
    scene_info = load_scene_info(scene_filename)

    # create new obstacle dataframe
    obs_df = scene_info.iloc[3:]

    # get number of obstacles
    num_obs = len(obs_df.index)
    
    # placeholders for object ids
    obs_id = []
    obs_name = []

    base_obs_path = 'assets/obstacles/'

    # iterate through obstacles
    for i in range(num_obs):
        # get obstacle info
        obs = list(obs_df.iloc[i])

        pos = obs[:3]
        scale = obs[3:-1]

        obs_filename = base_obs_path + Path(scene_filename).stem + '_obstacle_{}.urdf'.format(i)

        # save obstacle
        save_obstacle_urdf(obs_filename = obs_filename,
                        mass = 1,
                        inertia = [1,0,0,1,0,1],
                        mesh = "cube.obj",
                        scale = scale,
                        color = "red",
                        transparency = 0.8)

        obs_name.append(obs_filename)
        
        # load obstacle
        obs_id.append(p.loadURDF(obs_name[-1], basePosition=pos, useFixedBase=True))

        # set collision group
        p.setCollisionFilterGroupMask(obs_id[-1],-1, collision_group, enable_group_collision)

    return obs_name, obs_id

def save_obstacle_urdf(obs_filename = '',
                        mass = 1,
                        inertia = [1,0,0,1,0,1],
                        mesh = "cube.obj",
                        scale = [1,1,1],
                        color = "red",
                        transparency = 0.8):
    
    # more colors here: https://gist.github.com/naoki-mizuno/5e63a13597d5c5fe817c1600b829595e
    palette = {'red': [1, 0, 0],
                'green': [0, 1, 0],
                'blue': [0, 0, 1],
                'white': [1, 1, 1],
                'black': [1, 1, 1],
                'orange': [0.97, 0.451, 0.0235]}

    # set color
    rgba = palette[color] + [transparency]
    rgba = "{0} {1} {2} {3}".format(*rgba)

    # set inertial properties
    mass = str(mass)
    inertia = 'ixx="{0}" ixy="{1}" ixz="{2}" iyy="{3}" iyz="{4}" izz="{5}"'.format(*inertia)
    scale = "{0} {1} {2}".format(*scale)
    
    # set outputname
    name = Path(obs_filename).stem
    output_name = obs_filename

    # write file
    with open(output_name, "w") as f:
        text = '''<?xml version="1.0" ?>
            <robot name="{name}.urdf">
                <link name="baseLink">
                    <contact>
                        <lateral_friction value="0.5"/>
                        <inertia_scaling value="3.0"/>
                    </contact>
                    <inertial>
                        <origin rpy="0 0 0" xyz="0 0 0"/>
                        <mass value="{mass}"/>
                        <inertia {inertia}/>
                    </inertial>
                    <visual>
                        <origin rpy="0 0 0" xyz="0 0 0"/>
                        <geometry>
                            <mesh filename="{mesh}" scale="{scale}"/>
                        </geometry>
                        <material name="{color}">
                        <color rgba="{rgba}"/>
                        </material>
                    </visual>
                    <collision>
                        <origin rpy="0 0 0" xyz="0 0 0"/>
                        <geometry>
                            <!-- <mesh filename="{mesh}" scale="{scale}"/> --> 
                            <box size="{scale}"/>
                        </geometry>
                    </collision>
                </link>
            </robot>'''.format(name=name, mass=mass, inertia=inertia, mesh=mesh, scale=scale, color=color, rgba=rgba)
        f.write(text)
    return 

def load_robot(robot_filename, config, collision_group=None, enable_group_collision=False, rgba=None):
    print(robot_filename)
    robot_id = p.loadURDF(robot_filename, useFixedBase=True)
    
    if rgba is not None:
        for i in range(p.getNumJoints(robot_id)):
            p.changeVisualShape(
                robot_id,
                i,
                textureUniqueId=-1,
                rgbaColor=rgba)

    if collision_group is not None:
        for i in range(p.getNumJoints(robot_id)):
            p.setCollisionFilterGroupMask(robot_id, i, collision_group, enable_group_collision)

    joint_indices = list(range(len(config)))
    for joint_idx in joint_indices:
        p.resetJointState(robot_id, joint_idx, config[joint_idx])
    # p.setJointMotorControlArray(bodyIndex=robot_id, jointIndices=joint_indices, controlMode=p.POSITION_CONTROL, targetPositions=config)
    return robot_id

if __name__ == '__main__':

    # setup pybullet
    clid = p.connect(p.SHARED_MEMORY)

    if (clid < 0):
        p.connect(p.GUI, )#options="--background_color_red=1 --background_color_blue=1 --background_color_green=1")

    p.resetSimulation()
    p.setPhysicsEngineParameter(enableConeFriction=0)
    p.setAdditionalSearchPath(pybullet_data.getDataPath())

    p.setGravity(0, 0, -9.81)
    useRealTimeSimulation = 1
    p.setRealTimeSimulation(useRealTimeSimulation)
   
   
    # setup experiment
    task = 'manipulation'
    experiment= 'kinova_waitr_scenarios'
    planner = 'armour'
    basename = 'tasks/{}/{}/{}/'.format(task, experiment, planner)
    
    scene_dir = basename + 'scenes/'
    scene_names = os.listdir(scene_dir)
    scene_names = [scene_dir + name for name in scene_names]
    scene_names.sort()

    traj_dir = basename + 'traj/'
    traj_names = os.listdir(traj_dir)
    traj_names = [tra_dir + name for name in traj_names]
    traj_names.sort()

    
    # robot name
    base_robot_path = '../../urdfs/kinova/'
    robot_basename = base_robot_path + 'Kinova_Grasp_w_Tray_Gripper'
    # robot_basename = base_robot_path + 'kinova_arm/kinova_without_gripper'

    # flags
    save_sim = False
    reset_sim = True

    for scene_filename in scene_names:
        # create scene
        names, ids = create_static_scene(scene_filename, robot_basename, use_initial=True, use_current=False, use_bounding=False, use_goal=False)

        time.sleep(3)

        # create_recorder
        static_recorder = create_recorder(names, ids)


        static_recorder.add_keyframe()
        static_recorder.save('data/pkl_files/{}/{}/waitr_{}.pkl'.format(experiment, planner, Path(scene_filename).stem))

        # reset sim

        # remove obstacles
