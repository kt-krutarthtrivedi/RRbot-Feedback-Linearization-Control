%
% RBE502 - Spring 2023 | Programming Assignment 3
% Author: Krutarth Trivedi | ktrivedi@wpi.edu
%
clear; 
clc; 
close all;

global K;

% define symbols
syms theta1 theta2 theta1_dot theta2_dot theta1_ddot theta2_ddot t1 t2 'real'
syms m1 r1 l1 I1 m2 r2 l2 I2 g 'real' 'positive'

%--- Generate cubic polynomial trajectories for both the joints ------%
syms t;
timeVec = [1; t; t^2; t^3];
timeVec_dot = diff(timeVec);
timeVec_ddot = diff(timeVec_dot);

t0 = 0; 
tf = 10;
theta1_t0 = pi;
theta1_tf = 0;
theta2_t0 = pi/2;
theta2_tf = 0;
theta1_dot_t0 = 0;
theta1_dot_tf = 0;
theta2_dot_t0 = 0;
theta2_dot_tf = 0;

timeMat = [1, t0, t0^2, t0^3;
           0, 1, 2*t0, 3*t0^2;
           1, tf, tf^2, tf^3;
           0, 1, 2*tf, 3*tf^2];

theta1ConfigVec = [theta1_t0; theta1_dot_t0; theta1_tf; theta1_dot_tf];
theta2ConfigVec = [theta2_t0; theta2_dot_t0; theta2_tf; theta2_dot_tf];

theta1CoefficientVec = pinv(timeMat)*theta1ConfigVec;
theta2CoefficientVec = pinv(timeMat)*theta2ConfigVec;

q1 = theta1CoefficientVec' * timeVec;
q2 = theta2CoefficientVec' * timeVec;

q1_dot = theta1CoefficientVec' * timeVec_dot;
q2_dot = theta2CoefficientVec' * timeVec_dot;

q1_ddot = theta1CoefficientVec' * timeVec_ddot;
q2_ddot = theta2CoefficientVec' * timeVec_ddot;

%--------- Manipulator Equation Form ---------%
[M, C, G] = MEF();

% ------- Feedback Linearization --------%

% As calculated in Programming Assignment 2
A = [0         0    1.0000         0;
     0         0         0    1.0000;
     12.5769  -11.9611         0         0;
    -16.9227   46.1565         0         0];

B = [0         0;
     0         0;
     1.7250   -4.4345;
    -4.4345   14.8902];

desiredEigenValues = [-5,-10,-5,-10];
K = place(A, B, desiredEigenValues)

%--------------- Virtual control input ----------------------%
V = -K*([theta1; theta2; theta1_dot; theta2_dot] - [q1; q2; q1_dot; q2_dot]) + [q1_ddot; q2_ddot]

%-------------- The overall control law ---------------%
Tau = M*V + C*[q1_dot; q2_dot] + G

% -------- Simulation and Plotting ---------%

% simulate the system for 10 sec for given control inputs using ODE45
T = 10;
y0 = [deg2rad(200), deg2rad(125), 0 ,0];
[t,y] = ode45(@ode_rrbot, [0,T], y0);

%reconstruct the desired trajectories
q1_reconstruct = [];
q2_reconstruct = [];
q1_dot_reconstruct = [];
q2_dot_reconstruct= [];
q1_ddot_reconstruct = [];
q2_ddot_reconstruct= [];

M_reconstruct = [];
C_reconstruct = [];
G_reconstruct = [];

for i = 1:size(t)
    q1_reconstruct(i) = double(subs(q1, t(i,:)));
    q2_reconstruct(i) = double(subs(q2, t(i,:)));
    q1_dot_reconstruct(i) = double(subs(q1_dot, t(i,:)));
    q2_dot_reconstruct(i) = double(subs(q2_dot, t(i,:)));
    q1_ddot_reconstruct(i) = double(subs(q1_ddot, t(i,:)));
    q2_ddot_reconstruct(i) = double(subs(q2_ddot, t(i,:)));

    M_reconstruct(:,:,i) = double(subs(M, [theta1, theta2], [y(i,1), y(i,2)]));
    C_reconstruct(:,:,i) = double(subs(C, [theta1, theta2, theta1_dot, theta2_dot], [y(i,1), y(i,2), y(i,3), y(i,4)]));
    G_reconstruct = [G_reconstruct, double(subs(G, [theta1, theta2], [y(i,1), y(i,2)]))];
end

% reconstruct the applied control input
V = -K*(y' - [q1_reconstruct; q2_reconstruct; q1_dot_reconstruct; q2_dot_reconstruct]) + [q1_ddot_reconstruct; q2_ddot_reconstruct];

t1 = [];
t2 = [];

for i = 1:size(t)
    t1 = [t1, (M_reconstruct(1,:,i)*[V(1,i);V(2,i)]) + (C_reconstruct(1,:,i)*[q1_dot_reconstruct(1,i); q2_dot_reconstruct(1,i)]) + G_reconstruct(1,i)];
    t2 = [t2, (M_reconstruct(2,:,i)*[V(1,i);V(2,i)]) + (C_reconstruct(2,:,i)*[q1_dot_reconstruct(1,i); q2_dot_reconstruct(1,i)]) + G_reconstruct(2,i)];
end
% t1 = Tau(1,:);
% t2 = Tau(2,:);

% plot the trajectories
figure('Name','Trajectories', 'NumberTitle','off');
subplot(3,2,1)
plot(t,(y(:,1)),'b');
title('theta1')
xlabel('T');
ylabel('rad');
hold on;
plot(t,q1_reconstruct,'r');

subplot(3,2,2)
plot(t,(y(:,2)),'b')
title('theta2')
xlabel('T');
ylabel('rad');
hold on;
plot(t,q2_reconstruct,'r');

subplot(3,2,3)
plot(t,(y(:,3)),'b')
title('theta1-dot')
xlabel('T');
ylabel('rad/s');
hold on;
plot(t,q1_dot_reconstruct,'r');

subplot(3,2,4);
plot(t,(y(:,4)),'b');
title('theta2-dot')
xlabel('T');
ylabel('rad/s');
hold on;
plot(t,q2_dot_reconstruct,'r');

subplot(3,2,5);
plot(t,t1,'b');
title('t1')
xlabel('s');
ylabel('Nm');

subplot(3,2,6);
plot(t,t2,'b');
title('t2')
xlabel('s');
ylabel('Nm');

fprintf('Eigenvalues are selected such that: \n\n');
fprintf('Torque1: %.3f < u1 < %.3f \n\n',min(t1),max(t1));
fprintf('Torque2: %.3f < u2 < %.3f \n\n', min(t2),max(t2));