% Robot Controls - Feedback Linearization Control for RRbot Manipulator
% Author: Krutarth Trivedi | ktrivedi@wpi.edu

function [M,C,G] = MEF()
    syms theta1 theta2 theta1_dot theta2_dot theta1_ddot theta2_ddot t1 t2 'real'
    syms m1 r1 l1 I1 m2 r2 l2 I2 g 'real' 'positive'

    % physical parameters of the robot
    m1_ = 1; r1_ = 0.45; l1_ = 1; I1_ = 0.084;
    m2_ = 1; r2_ = 0.45; l2_ = 1; I2_ = 0.084;
    g_ = 9.81;

    % Symbolic EOM derived in programming assignment 1
    eom_1 = I1*theta1_ddot - t1 + (I2*(2*theta1_ddot + 2*theta2_ddot))/2 - (m2*(2*(l1*cos(theta1) + r2*cos(theta1 + theta2))*(r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot)^2 + l1*sin(theta1)*theta1_dot^2 - l1*theta1_ddot*cos(theta1) - r2*cos(theta1 + theta2)*(theta1_ddot + theta2_ddot)) - 2*(l1*sin(theta1) + r2*sin(theta1 + theta2))*(l1*cos(theta1)*theta1_dot^2 + l1*theta1_ddot*sin(theta1) + r2*sin(theta1 + theta2)*(theta1_ddot + theta2_ddot) + r2*cos(theta1 + theta2)*(theta1_dot + theta2_dot)^2)))/2 + (m1*(2*r1^2*theta1_ddot*cos(theta1)^2 + 2*r1^2*theta1_ddot*sin(theta1)^2))/2 - g*m2*(l1*sin(theta1) + r2*sin(theta1 + theta2)) - g*m1*r1*sin(theta1);
    eom_2 = (m2*(2*r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot)*(r2*cos(theta1 + theta2)*(theta1_dot + theta2_dot) + l1*cos(theta1)*theta1_dot) - 2*r2*cos(theta1 + theta2)*(l1*sin(theta1)*theta1_dot + r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot))*(theta1_dot + theta2_dot)))/2 - (m2*(2*r2*cos(theta1 + theta2)*(r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot)^2 + l1*sin(theta1)*theta1_dot^2 - l1*theta1_ddot*cos(theta1) - r2*cos(theta1 + theta2)*(theta1_ddot + theta2_ddot)) - 2*r2*sin(theta1 + theta2)*(l1*cos(theta1)*theta1_dot^2 + l1*theta1_ddot*sin(theta1) + r2*sin(theta1 + theta2)*(theta1_ddot + theta2_ddot) + r2*cos(theta1 + theta2)*(theta1_dot + theta2_dot)^2) + 2*r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot)*(r2*cos(theta1 + theta2)*(theta1_dot + theta2_dot) + l1*cos(theta1)*theta1_dot) - 2*r2*cos(theta1 + theta2)*(l1*sin(theta1)*theta1_dot + r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot))*(theta1_dot + theta2_dot)))/2 - t2 + (I2*(2*theta1_ddot + 2*theta2_ddot))/2 - g*m2*r2*sin(theta1 + theta2);
    
    % Calculate gravity vector
    G1 = subs(eom_1,[theta1_dot, theta2_dot, theta1_ddot, theta2_ddot, t1, t2],[0,0,0,0,0,0]);
    G2 = subs(eom_2,[theta1_dot, theta2_dot, theta1_ddot, theta2_ddot, t1, t2],[0,0,0,0,0,0]);
    G = [G1; G2];
    
    % Calculate Mass Matrix
    M11 = simplify(subs(eom_1 - G1,[theta1_dot, theta1_ddot, theta2_dot, theta2_ddot, t1, t2], [0,1,0,0,0,0]));
    M12 = simplify(subs(eom_1 - G1,[theta1_dot, theta1_ddot, theta2_dot, theta2_ddot, t1, t2], [0,0,0,1,0,0]));
    M21 = simplify(subs(eom_2 - G2,[theta1_dot, theta1_ddot, theta2_dot, theta2_ddot, t1, t2], [0,1,0,0,0,0]));
    M22 = simplify(subs(eom_2 - G2,[theta1_dot, theta1_ddot, theta2_dot, theta2_ddot, t1, t2], [0,0,0,1,0,0]));
    
    M = [M11, M12;
         M21, M22];
    
    % Calculate Coriolis Matrix
    C11 = simplify(subs(eom_1 - G1 - M11,[theta1_dot, theta1_ddot, theta2_dot, theta2_ddot, t1, t2], [1,0,0,0,0,0]));
    C12 = simplify(subs(eom_1 - G1 - M12,[theta1_dot, theta1_ddot, theta2_dot, theta2_ddot, t1, t2], [0,0,1,0,0,0]));
    C21 = simplify(subs(eom_2 - G2 - M21,[theta1_dot, theta1_ddot, theta2_dot, theta2_ddot, t1, t2], [1,0,0,0,0,0]));
    C22 = simplify(subs(eom_2 - G2 - M22,[theta1_dot, theta1_ddot, theta2_dot, theta2_ddot, t1, t2], [0,0,1,0,0,0]));
    
    C = [C11, C12;
        C21, C22];

    M = subs(M, [m1, r1, l1, I1, m2, r2, l2, I2, g], [m1_, r1_, l1_, I1_, m2_, r2_, l2_, I2_, g_]);
    C = subs(C, [m1, r1, l1, I1, m2, r2, l2, I2, g], [m1_, r1_, l1_, I1_, m2_, r2_, l2_, I2_, g_]);
    G = subs(G, [m1, r1, l1, I1, m2, r2, l2, I2, g], [m1_, r1_, l1_, I1_, m2_, r2_, l2_, I2_, g_]);
end