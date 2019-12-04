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

% clean-up
clear all, close all, clc;

% problem size
m = 1000;
n = 1000;

% randomize operator and offset
K = randi(3, m, n) - 2;
g = randi(11, n, 1) - 6;

% initialize starting point
x0 = 10 * randn(n, 1);
y0 = 10 * randn(m, 1);

% define primal/dual step sizes with    tau * sigma * L^2 < 1
[~, D] = svd(K, 'econ');
L = D(1, 1);
tau = sqrt((1 - 1e-4) / L ^ 2);
sigma = tau;

% function handles for F and G
lambda = 10;
F = @(y, c_flag) test_F(y, sigma, c_flag);
G = @(x, c_flag) test_G(x, g, lambda, tau, c_flag);

% perform optimization
theta = 1;
[x_star, y_star, primal_history, dual_history] = ...
    chambolle_pock(F, G, K, x0, y0, theta, tau, sigma);

% plot results
figure;
plot(1 : numel(primal_history(:, 1)), primal_history(:, 1), ...
    1 : numel(dual_history(:, 1)), dual_history(:, 1));
grid on;    axis tight;
legend('primal energy', 'dual energy');     xlabel('#iter');

%% F and G as local functions inside script

function [res1, res2, res3] = test_F(y, sigma, conjugate_flag)

res2 = 0;
if ~conjugate_flag
    res1 = 0.5 * sum(y .^ 2);
    res3 = 1 / (1 + sigma) * y;
else
    res1 = 0.5 * sum(y .^ 2);
    res3 = 1 / (1 + sigma) * y;
end

end

function [res1, res2, res3] = test_G(x, g, lambda, tau, conjugate_flag)

res2 = 0;
if ~conjugate_flag
    res1 = 0.5 * lambda * sum((x(:) - g(:)) .^ 2);
    res3 = 1 / (1 + (lambda * tau)) * (x + (lambda * tau) * g);
else
    res1 = lambda * (0.5 * sum((x / lambda) .^ 2) + ((x / lambda)' * g));
    res3 = [];
end

end