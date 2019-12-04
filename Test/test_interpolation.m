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

%% initialization
clear all, close all, clc;

% generate random image data
m = 3;
n = 4;
h = [1, 1];
img = randi(256, m, n) - 1;

%% display image data
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
colormap gray(256);
image(...
    'YData', [h(1) * (1/2), h(1) * (m - (1/2))], ...
    'XData', [h(2) * (1/2), h(2) * (n - (1/2))], ...
    'CData', img);
axis image;     set(gca, 'YDir', 'reverse');
colorbar;
xlabel('---y-->');      ylabel('<--x---');

% plot cell centered grid
omega = [0, m * h(1), 0, n * h(2)];
[x, y] = cell_centered_grid(omega, [m, n]);
p = [x(:), y(:)];
plot_grid(reshape(p, [m, n, 2]), 1, 'g--o');
title('image img with cell centered grid');

%% interpolate img over cell centered grid (... should yield img)
img_interpol = bilinear_interpolation(img, h, p);
img_interpol = reshape(img_interpol, [m, n]);
fprintf('||img - img_interpol|| = %.3e\n', norm(img(:) - img_interpol(:)));

%% randomize small displacement to test evaluate_displacement
u = 0.1 * randn(m * n, 2);
p_displaced = p + u;
img_u = evaluate_displacement(img, h, u);

% display results
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
colormap gray(256);

subplot(1, 2, 1);
image(...
    'YData', [h(1) * (1/2), h(1) * (m - (1/2))], ...
    'XData', [h(2) * (1/2), h(2) * (n - (1/2))], ...
    'CData', img);
axis image;     set(gca, 'YDir', 'reverse');
colorbar;
xlabel('---y-->');      ylabel('<--x---');
plot_grid(reshape(p_displaced, [m, n, 2]));
plot_grid(reshape(p, [m, n, 2]), 1, 'g--o');
hold on; quiver(p(:, 2), p(:, 1), u(:, 2), u(:, 1), 0, 'b'); hold off;
title('image img with displaced grid');

subplot(1, 2, 2);
image(...
    'YData', [h(1) * (1/2), h(1) * (m - (1/2))], ...
    'XData', [h(2) * (1/2), h(2) * (n - (1/2))], ...
    'CData', img_u);
axis image;     set(gca, 'YDir', 'reverse');
colorbar;
xlabel('---y-->');      ylabel('<--x---');
title('image img interpolated at displaced grid');