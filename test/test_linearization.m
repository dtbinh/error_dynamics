% Generating a test case for checking F, Q, H and R.

clear all;
close all;
clc;
addpath('..');

% Arbitrary state
p = [1;2;3];
v = [-2;3;1];
% q = [w, x, y, z]
% q = rpy2quat(10*pi/180, -5*pi/180, 230*pi/180);
q = [-0.424054855644904; 
0.002583605926050; 
0.097278948828693; 
0.900393031125399];
w = [0.2;-0.4;0.6];
kThrust = 8.55e-6;
kMoment = 0.016;
kInertia = [0.034756;0.045893;0.097700];
x = [p;v;q;w;kThrust;kMoment;kInertia];

u = [450; 500; 600; 324; 543; 800];

kMass = 1.6;
kArmLength = 0.2156;
kRotorPlaneOffset = 0.05;

dt = 0.01;

% Noise (here in STD, in C++ it's in COV)
% Process
process_noise_thrust = [2e-4; 2e-4; 1e-2];
process_noise_acceleration = 0.0 * ones(3,1);
process_noise_moment = 2e-3 * ones(3,1);

sigma_process_model = [repmat(process_noise_thrust, 6, 1);
    process_noise_acceleration;
    process_noise_moment];

% Measurement
meas_noise_pos = 0.005 * ones(3,1); 
meas_noise_att = 0.01 * ones(3,1); 

sigma_meas_model = [meas_noise_pos; meas_noise_att];


% Compute linear system matrices
F = compute_F_d(x, u, kMass, kArmLength, kRotorPlaneOffset, dt);
Q = compute_Q_d(x, u, kMass, kArmLength, kRotorPlaneOffset, sigma_process_model, dt);
H = compute_H_d(x, u, kMass, kArmLength, kRotorPlaneOffset);
R = compute_R_d(sigma_meas_model);

% Save matrices.
dlmwrite('F.txt', F(:)', 'precision', 16);
dlmwrite('Q.txt', Q(:)', 'precision', 16);
dlmwrite('H.txt', H(:)', 'precision', 16);
dlmwrite('R.txt', R(:)', 'precision', 16);