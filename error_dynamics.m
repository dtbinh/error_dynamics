% Calculate the linearized error dynamics.
clear all;
close all;
clc;

% General
syms dt real

% Firefly
syms kMass real
kNumRotors = 6;
syms kArmLength real;
kGravity = 9.807;
syms kRotorPlaneOffset real;
direction = [1,-1,1,-1,1,-1]';
motor_location = diag([kArmLength;kArmLength;1]) * ...
    [cos(30/180*pi), cos(90/180*pi), cos(150/180*pi), cos(210/180*pi), cos(270/180*pi), cos(330/180*pi);
    sin(30/180*pi), sin(90/180*pi), sin(150/180*pi), sin(210/180*pi), sin(270/180*pi), sin(330/180*pi);
    kRotorPlaneOffset, kRotorPlaneOffset, kRotorPlaneOffset, kRotorPlaneOffset, kRotorPlaneOffset, kRotorPlaneOffset];

% states / error states
p = sym_mat('p', 3, 1);
dp = sym_mat('dp', 3, 1);

v = sym_mat('v', 3, 1);
dv = sym_mat('dv', 3, 1);

q = sym_mat('q', 4, 1);
C = quat2rot(q);
dth = sym_mat('dth', 3, 1);
dC = eye(3) + skew(dth);

w = sym_mat('w', 3, 1);
dw = sym_mat('dw', 3, 1);

% parameters
kThrust = sym_mat('kThrust', 1, 1);
dkThrust = sym_mat('dkThrust', 1, 1);

kMoment = sym_mat('kMoment', 1, 1);
dkMoment = sym_mat('dkMoment', 1, 1);

kInertia = sym_mat('kInertia', 3, 1);
dkInertia = sym_mat('dkInertia', 3, 1);

kInertia_mat = diag(kInertia);
dkInertia_mat = diag(dkInertia);

% inputs
n_u = sym_mat('n_u', kNumRotors, 1);

% noise
% process
n_thrust = sym_mat('n_thrust', 3, kNumRotors);
sigma_n_thrust = repmat(sym_mat('sigma_n_thrust', 3, 1), 1, kNumRotors);

n_acc = sym_mat('n_acc', 3, 1);
sigma_n_acc = sym_mat('sigma_n_acc', 3, 1);

n_moment = sym_mat('n_moment', 3, 1);
sigma_n_moment = sym_mat('sigma_n_moment', 3, 1);

% measurements
n_p_m = sym_mat('n_p_m', 3, 1);
syms sigma_n_p_m real;

n_q_m = sym_mat('n_q_m', 3, 1);
syms sigma_n_q_m real;

% system states, input vector, noise and measurement vector

% system states, error states
x = [p; v; q; w; kThrust; kMoment; kInertia];
dx = [dp; dv; dth; dw; dkThrust; dkMoment; dkInertia];

% input vector
u = n_u;

% noise vector
n_process_model = [n_thrust(:); n_acc; n_moment];
sigma_process_model = [sigma_n_thrust(:); sigma_n_acc; sigma_n_moment];

n_meas_model = [n_p_m; n_q_m];
sigma_meas_model = [sigma_n_p_m * ones(3,1); sigma_n_q_m * ones(3,1)];

% Q_c, R_c
Q_c = diag(sigma_process_model.^2);
R_c = diag(sigma_meas_model.^2);

% state equations
thrusts = kThrust * [zeros(2,6); n_u.^2'];
thrusts_real = (kThrust + dkThrust) * [zeros(2,6); n_u.^2']  + n_thrust;

moments = - kMoment * repmat(direction', 3, 1) .* thrusts + ...
    cross(thrusts, motor_location);
moments_real = - (kMoment + dkMoment)* repmat(direction', 3, 1) .* thrusts_real ...
    + cross(thrusts_real, motor_location);

acc_pred = 1 / kMass * sum(thrusts,2) - C' * [0; 0; kGravity] - cross(w,v);
acc_real = 1 / kMass * sum(thrusts_real,2) - (C * dC)' * [0; 0; kGravity] ... 
    - cross(w + dw, v + dv) + n_acc;

w_dot_pred = inv(kInertia_mat) * (sum(moments,2) - ...
    cross(w, kInertia_mat * w));
w_dot_real = inv(kInertia_mat + dkInertia_mat) * (sum(moments_real,2) - ...
    cross(w + dw, (kInertia_mat + dkInertia_mat) * (w + dw))) + n_moment;

Omega = [0, -w';
    w, -skew(w)];

% nominal continuous time model
f_state = [
    C * v;
    acc_pred;
    0.5 * Omega * q;
    w_dot_pred;
    0;                          % kThrust
    0;                          % kMoment
    zeros(3,1);                 % kInertia
    ];
f_state = simplify(f_state);

% error states
dp_dot = (C * dC) * (v + dv) - C * v;
dv_dot = acc_real - acc_pred;
dth_dot = -skew(w) * dth + dw;
dw_dot = w_dot_real - w_dot_pred;
dkThrust_dot = 0;
dkMoment_dot = 0;
dkInertia_dot = zeros(3,1);

dx_dot = [dp_dot; dv_dot; dth_dot; dw_dot; dkThrust_dot; dkMoment_dot; dkInertia_dot];
dx_dot = simplify(dx_dot);

% linearized error system matrix F_c
F_c = jacobian(dx_dot, dx);
% remove noise (zero mean) and higher order error terms
F_c = subs(F_c, dx, zeros(size(dx)));
F_c = subs(F_c, n_process_model, zeros(size(n_process_model)));
F_c = simplify(F_c);
% discrete version
F_d = eye(size(F_c, 1)) + F_c * dt; % zero order hold
F_d = simplify(F_d);
disp('calculated Fd');


% linearized noise matrix
G_c = jacobian(dx_dot, n_process_model);
% remove remaing error states (zero steady state) and higher order noise terms
G_c = subs(G_c, dx, zeros(size(dx)));
G_c = subs(G_c, n_process_model, zeros(size(n_process_model)));
G_c = simplify(G_c);
% discrete version
Q_d = dt * G_c * Q_c * G_c';
Q_d = simplify(Q_d);
disp('calculated Qd');

% calculate update equations
dz_p = dp + n_p_m;
dz_th = dth + n_q_m;

dz = [dz_p; dz_th];

% linearized output matrix
H = jacobian(dz, dx);
% remove noise (zero mean) and higher order error terms
H = subs(H, dx, zeros(size(dx)));
H = subs(H, n_meas_model, zeros(size(n_meas_model)));
H = simplify(H);
% discrete version
H_d = H;
disp('calculated Hd');

% measurement noise transition matrix
V = jacobian(dz, n_meas_model);
% remove noise (zero mean) and higher order error terms
V = subs(V, dx, zeros(size(dx)));
V = subs(V, n_meas_model, zeros(size(n_meas_model)));
V = simplify(V);
disp('calculated Vd');

% measurement noise
R_d = R_c;

% Generate matlab functions.
disp('generating code to compute Fd');
matlabFunction(F_d, 'file', 'compute_F_d', 'vars', {x, n_u, kMass, kArmLength, kRotorPlaneOffset, dt});
disp('generating code to compute Qd');
matlabFunction(Q_d, 'file', 'compute_Q_d', 'vars', {x, n_u, kMass, kArmLength, kRotorPlaneOffset, sigma_process_model, dt});
disp('generating code to compute Hd');
matlabFunction(H_d, 'file', 'compute_H_d', 'vars', {x, n_u, kMass, kArmLength, kRotorPlaneOffset});
disp('generating code to compute Rd');
matlabFunction(R_d, 'file', 'compute_R_d', 'vars', {sigma_meas_model});