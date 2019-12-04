%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------

% add project to path
addpath(genpath('..'));

%% 3D (FOR TESTING PROJECTION)

clear all, close all, clc;

[X, Y, Z] = ndgrid(-1 : 0.02 : 1, -1 : 0.02 : 1, -1 : 0.02 : 1);
P = [X(:), Y(:), Z(:)];
idx = abs(sum(abs(P), 2) - 1) <= 1e-3;
P = P(idx, :);

figure;
scatter3(P(:, 1), P(:, 2), P(:, 3), 'r.');
axis equal;

Q = [...
    +1,  0,  0; ...
    -1,  0,  0; ...
    0, +1,  0; ...
    0, -1,  0; ...
    0,  0, +1; ...
    0,  0, -1; ...
    ];

L = [...
    1, 3, 5; ...
    1, 3, 6; ...
    1, 4, 5; ...
    1, 4, 6; ...
    2, 3, 5; ...
    2, 3, 6; ...
    2, 4, 5; ...
    2, 4, 6; ...
    ];

q = rand(3, 1);
while sum(abs(q)) <= 1
    q = rand(3, 1);
end
r = l1ball_projection(q);
fprintf('Distance ||q - Proj(q)||_2 = %.4f\n', ...
    sqrt(sum((q - r) .^ 2)));

% compare to minimum distance between q and a grid of surface points P
surf_dist = sqrt(sum((P - q') .^ 2, 2));
fprintf('Minimum distance to surface grid points = %.4f\n', ...
    min(surf_dist));

hold on;
for i = 1 : size(L, 1)
    tmp = Q([L(i, :), L(i, 1)], :);
    plot3(tmp(:, 1), tmp(:, 2), tmp(:, 3), 'k', 'LineWidth', 2);
end
plot3([q(1), r(1)], [q(2), r(2)], [q(3), r(3)], 'b-o');
hold off;

%% ND (FOR TESTING TIME COMPLEXITY)

clear all, close all, clc;

N = 5 * 10 .^ (0 : 5);
T = zeros(size(N));

for n = 1 : numel(N)
    for i = 1 : 10
        q = randn(1, N(n));
        q = 1.1 * q / sum(abs(q));
        tic;
        l1ball_projection(q);
        T(n) = T(n) + toc;
    end
    T(n) = T(n) / 10;
end

figure;
loglog(N, T, '-o');
grid on;