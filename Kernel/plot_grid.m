function plot_grid(g, s, style, f)
%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------
% IN:
%   g       ~ m x n x 2             grid to plot
%   s       ~ 1 x 1                 spacing for grid plotting
%   style   ~ string                string for plotting style
%   f       ~ figure handle         figure to plot in
%--------------------------------------------------------------------------

% set active figure
if nargin < 4
    f = gcf;
end
if nargin < 3
    style = 'm-o';
end
if nargin < 2
    s = 1;
end
figure(f);
hold on;

% fetch grid size
[m, n, ~] = size(g);

% exchange x-/y-components in grid g in order to use plot(..) later on
%   ~> plot(..) always interpretes the horizontal as the first axis!
g = cat(3, g(:, :, 2), g(:, :, 1));

%
IDX1 = round(linspace(1, m, m / s));
IDX2 = round(linspace(1, n, n / s));

% plot horizontal grid lines
for i = IDX1
    plot(g(i, :, 1), g(i, :, 2), style);
end

% plot vertical grid lines
for j = IDX2
    plot(g(:, j, 1), g(:, j, 2), style);
end

hold off;

end