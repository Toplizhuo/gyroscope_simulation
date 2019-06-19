% Script written by:
% Zhuo Li (zhuol7@student.unimelb.edu.au)
% The University of Melbourne

%% Define the seperate rotation matrix

syms a(t) b(t) c(t) d(t);  

R01 = [cos(a(t)) -sin(a(t)) 0; sin(a(t)) cos(a(t)) 0; 0 0 1];
R10 = [cos(a(t)) sin(a(t)) 0; -sin(a(t)) cos(a(t)) 0; 0 0 1];
R12 = [1 0 0; 0 cos(b(t)) -sin(b(t)); 0 sin(b(t)) cos(b(t))];
R21 = [1 0 0; 0 cos(b(t)) sin(b(t)); 0 -sin(b(t)) cos(b(t))];
R23 = [cos(c(t)) -sin(c(t)) 0; sin(c(t)) cos(c(t)) 0; 0 0 1];
R32 = [cos(c(t)) sin(c(t)) 0; -sin(c(t)) cos(c(t)) 0; 0 0 1];
R34 = [cos(d(t)) -sin(d(t)) 0; sin(d(t)) cos(d(t)) 0; 0 0 1];
R43 = [cos(d(t)) sin(d(t)) 0; -sin(d(t)) cos(d(t)) 0; 0 0 1];

w_1_1 = [0; 0; diff(a(t), t)];
w_21_2 = [diff(b(t), t); 0; 0];
w_32_3 = [0; 0; diff(c(t), t)];
w_43_4 = [0; 0; diff(d(t), t)];

w_2_2 = w_21_2 + R21 * w_1_1;
w_3_3 = R32 * w_2_2 + w_32_3;
w_r_3 = w_3_3 + w_43_4;
w_4_4 = w_43_4 + R43 * w_3_3;

%% Define Variables

syms R_f R_r R_i real      % R_f, R_r, R_i is the radius of the frame, rotor and the frame ring
syms M_f M_r real          % M_f, M_r is the mass of the frame and rotor
syms H_t H_r real          % H_t, H_r is the total height and the height of the rotor

R_f = 32.75; % radius of the frame
R_r = 28.50; % ratius of the rotor
R_i = 3.25; % ratius of the frame ring
M_f = 24;   % mass of the frame
M_r = 45;   % mass of the rotor
H_t = 94;   % total height
H_r = 65;   % height of the rotor
g = 9.80;      % acceleration of gravity
L_OG=47;    % distance between point O and G

%% Calculate the inertia tensor of frame on point G & O

% the mass of pole
M_pole=((pi*R_i^2)*(H_t-2*R_f)*M_f)/(2*(pi*R_i^2)*(pi*2*R_f)+(pi*R_i^2)*(H_t-2*R_f));
% the mass of ring
M_ring=((pi*R_i^2)*(pi*2*R_f)*M_f)/(2*(pi*R_i^2)*(pi*2*R_f)+(pi*R_i^2)*(H_t-2*R_f));

% the inertia tensor of rotor
I_G_r=M_r.*[(3*R_r^2+H_r^2)/12 0 0; 0 (3*R_r^2+H_r^2)/12 0; 0 0 (R_r^2)/2];

% the inertia tensor of the Part a on point G
I_G_a = (M_pole/12).*[(3*R_i^2+H_t^2) 0 0; 0 (3*R_i^2+H_t^2) 0; 0 0 (R_i^2)*6];
% the inertia tensor of the Part b on point G
I_G_b = M_ring.*[(5/8)*R_i^2+(1/2)*R_f^2 0 0; 0 (5/8)*R_i^2+(1/2)*R_f^2 0; 0 0 (3/4)*R_i^2+R_f^2];
% the inertia tensor of the Part c on point G
I_G_c = M_ring.*[(5/8)*R_i^2+(1/2)*R_f^2 0 0; 0 (3/4)*R_i^2+R_f^2 0; 0 0 (5/8)*R_i^2+(1/2)*R_f^2];
% the inertia tensor of frame on point G
I_G_f = I_G_a + I_G_b + I_G_c;

% the inertia tensor of frame on point O
I_O_f = I_G_f + M_f.*[L_OG^2 0 0; 0 L_OG^2 0; 0 0 0]; 

%% Newton-Euler equations

syms F_G_3_x F_G_3_y F_G_3_z M_G_3_x M_G_3_y F_O_3_x F_O_3_y F_O_3_z;

F_G_3 = [F_G_3_x; F_G_3_y; F_G_3_z];
M_G_3 = [M_G_3_x; M_G_3_y; 0];
F_O_3 = [F_O_3_x; F_O_3_y; F_O_3_z];

% Newton-Euler equations for rotor - Linear

r_OG_3 = [0; 0; L_OG];
r_OG_3_1 = diff(r_OG_3, t) + cross(w_3_3, r_OG_3);
r_OG_3_2 = diff(r_OG_3_1, t) + cross(w_3_3, r_OG_3_1);

G_r_0 = [0; 0; -M_r*g];
G_r_3 = R32*R21*R10*G_r_0;
eqns_r_L = simplify(F_G_3 + G_r_3) == simplify(M_r.*r_OG_3_2);
F_G_3_x = solve(eqns_r_L(1), F_G_3_x);
F_G_3_y = solve(eqns_r_L(2), F_G_3_y);
F_G_3_z = solve(eqns_r_L(3), F_G_3_z);

F_G_3 = [F_G_3_x; F_G_3_y; F_G_3_z]; 	

% Newton-Euler equations for rotor - Angular

h_r_4_G = I_G_r * w_4_4;
h_r_3_G = R34 * h_r_4_G;
eqns_r_A = simplify(M_G_3) == simplify(diff(h_r_3_G, t) + cross(w_3_3, h_r_3_G));
M_G_3_x = solve(eqns_r_A(1), M_G_3_x);
M_G_3_y = solve(eqns_r_A(2), M_G_3_y);
M_G_3 = [M_G_3_x; M_G_3_y; 0];    	                       
                         
% Newton-Euler equations for frame - Linear
                       
G_f_0 = [0; 0; -M_f * g];
G_f_3 = R32 * R21 * R10 * G_f_0;
EQ_f_L = simplify(-F_G_3 + F_O_3 + G_f_3) == simplify(M_f .* r_OG_3_2);
F_O_3_x = solve(EQ_f_L(1), F_O_3_x);
F_O_3_y = solve(EQ_f_L(2), F_O_3_y);
F_O_3_z = solve(EQ_f_L(3), F_O_3_z);
F_O_3 = [F_O_3_x; F_O_3_y; F_O_3_z];

% Newton-Euler equations for frame - Angular   

h_f_3_O = I_O_f * w_3_3;
equs_f_A = simplify(-M_G_3 + cross(r_OG_3, -F_G_3)) == simplify(diff(h_f_3_O, t) + cross(w_3_3, h_f_3_O));
eqns = simplify([eqns_r_A(3); equs_f_A]);

% Replaces symfun-type variables with sym-type variables

syms a_0 b_0 c_0 d_0 real;           % a, b, c, d equal to alpha, beta, gamma, delta in the report 
syms a_1 b_1 c_1 d_1 real;   % x_1 is the first time derivatives of x 
syms a_2 b_2 c_2 d_2 real;   % x_2 is the Second time derivatives of x

eqns = subs(eqns, [diff(a(t), t, t), diff(b(t), t, t), diff(c(t), t, t), diff(d(t), t, t)], [a_2, b_2, c_2, d_2]);
eqns = subs(eqns, [diff(a(t), t), diff(b(t), t), diff(c(t), t), diff(d(t), t)], [a_1, b_1, c_1, d_1]); 
eqns = subs(eqns, [a(t), b(t), c(t), d(t)], [a_0, b_0, c_0, d_0]);

[A, B] = equationsToMatrix(eqns, [a_2, b_2, c_2, d_2]);

expression = simplify(A\B)