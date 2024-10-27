# EEE4022S-Project-2024-NXSMPI001
This GitHub is repo contains the simulink files and generated C++ ROS nodes for my EEE4022S final year project. This project investigated the application of Nonlinear Model Predictive Control for trajectory tracking in autonomous surface vessels. To achieve this, an NMPC controller was developed for deployment on the Subsea Tech Catarob. 

## Simulink Implementation
The NMPC controller was developed and simulated in MATLAB simulink. To simulate the closed loop system, go to the [Folder Name](Simulink/) folder and ...
1.  Open and run the 'ControlParams.mlx' script in MATLAB
2.  Open and run the 'nmpcSim.slx' model in Simulink
The final controller for implementation in ROS can be found in the Simulink models 'nmpcRect.slx' and 'nmpcSine.slx'. The two ROS nodes for rectangular and sinusodial reference paths, respectively, are generated from these models using the ROS toolbox and Simulink Coder.

## ROS Implementation
The two generated C++ nodes can be found in the ROS2 workspace under the [Folder Name](Ros2_ws/) folder. To run the, for example, the '/nmpcRect' node ...
1.  Copy the  [Folder Name](Ros2_ws/src/) folder to your own ROS2 workspace.
2.  Run `colcon build` to build the workspace.
3.  Run `source <path_to_ws>/install/setup.bash` to source the workspace.
4.  Run `<path_to_ws>/install/nmpcrect/lib/nmpcrect/nmpcRect` to start the ROS node.


