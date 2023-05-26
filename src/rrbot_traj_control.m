% Robot Controls - Feedback Linearization Control for RRbot Manipulator
% Author: Krutarth Trivedi | ktrivedi@wpi.edu

clear; close; clc;

% ROS Setup
rosinit;

% global variable to store Gain Matrix
global K;
% K = 1.0e+03 * [1.1604   -0.2123    0.1151    0.0120;
%     0.3281   -0.0221    0.0336    0.0070];


j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');
JointStates = rossubscriber('/rrbot/joint_states');

tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);

tau1.Data = 0;
tau2.Data = 0;

send(j1_effort,tau1);
send(j2_effort,tau2);

client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(200), deg2rad(125)];
resp = call(client,req,'Timeout',3);

tic;
t = 0;

sample_theta1 = [];
sample_theta2 = [];
sample_theta1_dot = [];
sample_theta2_dot = [];
sample_time = [];
sample_tau1 = [];
sample_tau2 = [];

while(t < 10)

    t = toc;
    % read the joint states
    jointData = receive(JointStates);

    % implementing the trajectory feedback controller  
    theta1 = jointData.Position(1,1);
    theta2 = jointData.Position(2,1);
    theta1_dot = jointData.Velocity(1,1);
    theta2_dot = jointData.Velocity(2,1);

    %Get the trajectory

    X = [theta1, theta2, theta1_dot, theta2_dot];

    %--- Generate cubic polynomial trajectories for both the joints ------%
     
    q1 = (pi*t^3)/500 - (3*pi*t^2)/100 - (6189958033024885*t)/10141204801825835211973625643008 + pi;
    q2 = (pi*t^3)/1000 - (3*pi*t^2)/200 - (6189958033024885*t)/20282409603651670423947251286016 + pi/2;
    q1_dot = (3*pi*t^2)/500 - (3*pi*t)/50 - 6189958033024885/10141204801825835211973625643008;
    q2_dot = (3*pi*t^2)/1000 - (3*pi*t)/100 - 6189958033024885/20282409603651670423947251286016;
    q1_ddot = (3*pi*t)/250 - (3*pi)/50; 
    q2_ddot = (3*pi*t)/500 - (3*pi)/100;
    
    %--------- Manipulator Equation Form ---------%
    M = [(9*cos(theta2))/10 + 1573/1000, (9*cos(theta2))/20 + 573/2000;
        (9*cos(theta2))/20 + 573/2000,                      573/2000];
     
     
    C =[- (9*cos(theta2))/10 - 1573/1000, - (9*cos(theta2))/20 - (9*sin(theta2))/20 - 573/2000;
       (9*sin(theta2))/20 - (9*cos(theta2))/20 - 573/2000, -573/2000];
     
     
    G =[- (8829*sin(theta1 + theta2))/2000 - (28449*sin(theta1))/2000;
                                -(8829*sin(theta1 + theta2))/2000];
    
    
    %--------------- Virtual control input ----------------------%
    V = -K*([theta1; theta2; theta1_dot; theta2_dot] - [q1; q2; q1_dot; q2_dot]) + [q1_ddot; q2_ddot];

    %-------------- The overall control law ---------------%
    Tau = M*V + C*[theta1_dot; theta2_dot] + G;
 
    tau1.Data = Tau(1);
    tau2.Data = Tau(2);
    
    send(j1_effort,tau1);
    send(j2_effort,tau2);

    % sample the time, joint state values, and calculated torques here to be plotted at the end
    sample_time = [sample_time, t];
    sample_theta1 = [sample_theta1, theta1];
    sample_theta2 = [sample_theta2, theta2];
    sample_theta1_dot = [sample_theta1_dot, theta1_dot];
    sample_theta2_dot = [sample_theta2_dot, theta2_dot];
    sample_tau1 = [sample_tau1, Tau(1)];
    sample_tau2 = [sample_tau2, Tau(2)];
    
end

tau1.Data = 0;
tau2.Data = 0;

send(j1_effort,tau1);
send(j2_effort,tau2);

% disconnect from roscore
rosshutdown;

% plot the trajectories
figure('Name','State Trajectories', 'NumberTitle','off');
subplot(3,2,1)
plot(sample_time,sample_theta1,'b');
title('theta1')
xlabel('T');
ylabel('rad');

subplot(3,2,2)
plot(sample_time,sample_theta2,'b')
title('theta2')
xlabel('T');
ylabel('rad');

subplot(3,2,3)
plot(sample_time,sample_theta1_dot,'b')
title('theta1-dot')
xlabel('T');
ylabel('rad/s');

subplot(3,2,4);
plot(sample_time,sample_theta2_dot,'b');
title('theta2-dot')
xlabel('T');
ylabel('rad/s');

subplot(3,2,5);
plot(sample_time,sample_tau1,'b');
title('tau1')
xlabel('s');
ylabel('Nm');

subplot(3,2,6);
plot(sample_time,sample_tau2,'b');
title('tau2')
xlabel('s');
ylabel('Nm');