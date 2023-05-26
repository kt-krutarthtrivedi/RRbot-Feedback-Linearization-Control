%
% RBE502 - Spring 2023 | Programming Assignment 3
% Author: Krutarth Trivedi | ktrivedi@wpi.edu
%
function dX = ode_rrbot(t,X)

global K;

% physical parameters of the robot
m1 = 1; r1 = 0.45; l1 = 1; I1 = 0.084;
m2 = 1; r2 = 0.45; l2 = 1; I2 = 0.084;
g = 9.81;

%--------------- Get the current state values ---------------%
dX = zeros(4,1);
X = num2cell(X);
[theta1, theta2, theta1_dot, theta2_dot] = deal(X{:});

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
Tau = M*V + C*[q1_dot; q2_dot] + G;

t1 = Tau(1);
t2 = Tau(2);

dX(1) = theta1_dot;
dX(2) = theta2_dot;
dX(3) = (I2*t1 - I2*t2 + m2*r2^2*t1*cos(theta1 + theta2)^2 - m2*r2^2*t2*cos(theta1 + theta2)^2 + m2*r2^2*t1*sin(theta1 + theta2)^2 - m2*r2^2*t2*sin(theta1 + theta2)^2 + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2^2*r2^3*cos(theta1)*sin(theta1 + theta2)^3*theta1_dot^2 - l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)^3*theta1_dot^2 + l1*m2^2*r2^3*cos(theta1)*sin(theta1 + theta2)^3*theta2_dot^2 - l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)^3*theta2_dot^2 - l1*m2*r2*t2*sin(theta1)*sin(theta1 + theta2) + g*l1*m2^2*r2^2*sin(theta1)*cos(theta1 + theta2)^2 - l1*m2*r2*t2*cos(theta1)*cos(theta1 + theta2) - l1^2*m2^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta1_dot^2 + I2*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta1_dot^2 - I2*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta1_dot^2 + I2*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta2_dot^2 - I2*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta2_dot^2 + g*m1*m2*r1*r2^2*sin(theta1)*cos(theta1 + theta2)^2 + g*m1*m2*r1*r2^2*sin(theta1)*sin(theta1 + theta2)^2 - g*l1*m2^2*r2^2*cos(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2) + 2*l1*m2^2*r2^3*cos(theta1)*sin(theta1 + theta2)^3*theta1_dot*theta2_dot - 2*l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)^3*theta1_dot*theta2_dot + l1*m2^2*r2^3*cos(theta1)*cos(theta1 + theta2)^2*sin(theta1 + theta2)*theta1_dot^2 + l1*m2^2*r2^3*cos(theta1)*cos(theta1 + theta2)^2*sin(theta1 + theta2)*theta2_dot^2 - l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2)^2*theta1_dot^2 - l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2)^2*theta2_dot^2 - l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*cos(theta1 + theta2)^2*theta1_dot^2 + l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*sin(theta1 + theta2)^2*theta1_dot^2 + l1^2*m2^2*r2^2*cos(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta1_dot^2 + 2*l1*m2^2*r2^3*cos(theta1)*cos(theta1 + theta2)^2*sin(theta1 + theta2)*theta1_dot*theta2_dot - 2*l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2)^2*theta1_dot*theta2_dot + 2*I2*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta1_dot*theta2_dot - 2*I2*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta1_dot*theta2_dot)/(I1*I2 + I2*l1^2*m2*cos(theta1)^2 + I2*m1*r1^2*cos(theta1)^2 + I2*l1^2*m2*sin(theta1)^2 + I2*m1*r1^2*sin(theta1)^2 + I1*m2*r2^2*cos(theta1 + theta2)^2 + I1*m2*r2^2*sin(theta1 + theta2)^2 + l1^2*m2^2*r2^2*cos(theta1)^2*sin(theta1 + theta2)^2 + l1^2*m2^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*cos(theta1)^2*cos(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*cos(theta1)^2*sin(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*sin(theta1)^2*sin(theta1 + theta2)^2 - 2*l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2));
dX(4) = (I1*t2 - I2*t1 + I2*t2 - m2*r2^2*t1*cos(theta1 + theta2)^2 + m2*r2^2*t2*cos(theta1 + theta2)^2 - m2*r2^2*t1*sin(theta1 + theta2)^2 + m2*r2^2*t2*sin(theta1 + theta2)^2 + l1^2*m2*t2*cos(theta1)^2 + m1*r1^2*t2*cos(theta1)^2 + l1^2*m2*t2*sin(theta1)^2 + m1*r1^2*t2*sin(theta1)^2 - I2*g*l1*m2*sin(theta1) - I2*g*m1*r1*sin(theta1) + I1*g*m2*r2*sin(theta1 + theta2) - l1*m2^2*r2^3*cos(theta1)*sin(theta1 + theta2)^3*theta1_dot^2 + l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)^3*theta1_dot^2 - l1^3*m2^2*r2*cos(theta1)^3*sin(theta1 + theta2)*theta1_dot^2 + l1^3*m2^2*r2*sin(theta1)^3*cos(theta1 + theta2)*theta1_dot^2 - l1*m2^2*r2^3*cos(theta1)*sin(theta1 + theta2)^3*theta2_dot^2 + l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)^3*theta2_dot^2 - l1*m2*r2*t1*sin(theta1)*sin(theta1 + theta2) + 2*l1*m2*r2*t2*sin(theta1)*sin(theta1 + theta2) - g*l1*m2^2*r2^2*sin(theta1)*cos(theta1 + theta2)^2 + g*l1^2*m2^2*r2*cos(theta1)^2*sin(theta1 + theta2) - l1*m2*r2*t1*cos(theta1)*cos(theta1 + theta2) + 2*l1*m2*r2*t2*cos(theta1)*cos(theta1 + theta2) + 2*l1^2*m2^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta1_dot^2 + l1^2*m2^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta2_dot^2 - g*l1^2*m2^2*r2*cos(theta1)*sin(theta1)*cos(theta1 + theta2) - I1*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta1_dot^2 + I1*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta1_dot^2 - I2*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta1_dot^2 + I2*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta1_dot^2 - I2*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta2_dot^2 + I2*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta2_dot^2 - g*m1*m2*r1*r2^2*sin(theta1)*cos(theta1 + theta2)^2 + g*m1*m2*r1^2*r2*cos(theta1)^2*sin(theta1 + theta2) - g*m1*m2*r1*r2^2*sin(theta1)*sin(theta1 + theta2)^2 + g*m1*m2*r1^2*r2*sin(theta1)^2*sin(theta1 + theta2) + l1^3*m2^2*r2*cos(theta1)^2*sin(theta1)*cos(theta1 + theta2)*theta1_dot^2 + g*l1*m2^2*r2^2*cos(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2) - l1^3*m2^2*r2*cos(theta1)*sin(theta1)^2*sin(theta1 + theta2)*theta1_dot^2 - 2*l1*m2^2*r2^3*cos(theta1)*sin(theta1 + theta2)^3*theta1_dot*theta2_dot + 2*l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)^3*theta1_dot*theta2_dot - l1*m2^2*r2^3*cos(theta1)*cos(theta1 + theta2)^2*sin(theta1 + theta2)*theta1_dot^2 - l1*m2^2*r2^3*cos(theta1)*cos(theta1 + theta2)^2*sin(theta1 + theta2)*theta2_dot^2 + l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2)^2*theta1_dot^2 + l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2)^2*theta2_dot^2 + 2*l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*cos(theta1 + theta2)^2*theta1_dot^2 + l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*cos(theta1 + theta2)^2*theta2_dot^2 - 2*l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*sin(theta1 + theta2)^2*theta1_dot^2 - l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*sin(theta1 + theta2)^2*theta2_dot^2 - 2*l1^2*m2^2*r2^2*cos(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta1_dot^2 - l1^2*m2^2*r2^2*cos(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta2_dot^2 - 2*l1*m2^2*r2^3*cos(theta1)*cos(theta1 + theta2)^2*sin(theta1 + theta2)*theta1_dot*theta2_dot + 2*l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2)^2*theta1_dot*theta2_dot + 2*l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*cos(theta1 + theta2)^2*theta1_dot*theta2_dot - 2*l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*sin(theta1 + theta2)^2*theta1_dot*theta2_dot - g*l1*m1*m2*r1*r2*sin(theta1)^2*sin(theta1 + theta2) - 2*l1^2*m2^2*r2^2*cos(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta1_dot*theta2_dot - l1*m1*m2*r1^2*r2*cos(theta1)^3*sin(theta1 + theta2)*theta1_dot^2 + l1*m1*m2*r1^2*r2*sin(theta1)^3*cos(theta1 + theta2)*theta1_dot^2 + 2*l1^2*m2^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta1_dot*theta2_dot - 2*I2*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta1_dot*theta2_dot + 2*I2*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta1_dot*theta2_dot - g*l1*m1*m2*r1*r2*cos(theta1)*sin(theta1)*cos(theta1 + theta2) + l1*m1*m2*r1^2*r2*cos(theta1)^2*sin(theta1)*cos(theta1 + theta2)*theta1_dot^2 - l1*m1*m2*r1^2*r2*cos(theta1)*sin(theta1)^2*sin(theta1 + theta2)*theta1_dot^2)/(I1*I2 + I2*l1^2*m2*cos(theta1)^2 + I2*m1*r1^2*cos(theta1)^2 + I2*l1^2*m2*sin(theta1)^2 + I2*m1*r1^2*sin(theta1)^2 + I1*m2*r2^2*cos(theta1 + theta2)^2 + I1*m2*r2^2*sin(theta1 + theta2)^2 + l1^2*m2^2*r2^2*cos(theta1)^2*sin(theta1 + theta2)^2 + l1^2*m2^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*cos(theta1)^2*cos(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*cos(theta1)^2*sin(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*sin(theta1)^2*sin(theta1 + theta2)^2 - 2*l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2));

end