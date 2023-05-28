# RRbot-Feedback-Linearization-Control

## Aim
The aim is to design a feedback linearization control law for the RRbot and perform trajectory tracking.

## Design Procedure
The desired cubic polynomial trajectories are calculated based on the state parameters given for the initial and final conditions.

The control law is designed in MATLAB, and the parameters are tuned to achieve desirable state trajectories and control input trajectories. Once the results are satisfactory, further experiments are performed in Gazebo.

## Results

- The implemented control law performed as expected to perform trajectory tracking within the given limited torque constraints. 
- The state and control input trajectories are captured as follows: 

<img width="681" alt="Screenshot 2023-05-27 at 11 21 23 PM" src="https://github.com/kt-krutarthtrivedi/RRbot-Feedback-Linearization-Control/assets/134632027/e9465564-f0d1-4571-a8ea-21bd93f899c9">



## Report

Please find a detailed description of the control law design and results in the [report](https://github.com/kt-krutarthtrivedi/RRbot-Feedback-Linearization-Control/blob/main/media/Report.pdf).
