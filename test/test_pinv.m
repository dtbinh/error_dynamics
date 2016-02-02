%% Groundtruth to test pseudo inverse implementation

% 1. Allocation matrix test
s = sin(30 * pi / 180);
c = cos(30 * pi / 180);
A =[ s, 1, s, -s, -1, -s;
    -c, 0, c, c, 0, -c;
    -1, 1, -1, 1, -1, 1;
     1, 1, 1, 1, 1, 1 ];
arm_length = 0.2156;
c_t = 8.55e-6;
c_m = 0.016;
K = diag([arm_length * c_t, arm_length * c_t, c_t * c_m, c_t]);

Alloc_inv = pinv(A) * inv(K);

dlmwrite('Alloc_inv.txt', Alloc_inv(:)', 'precision', 16);
