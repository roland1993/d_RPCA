function [img_u, dimg_u] = evaluate_displacement(img, h, u, omega, s_grid)
%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------
% IN:
%   img     ~ m x n             image
%   h       ~ 2 x 1             image step sizes [h_x, h_y]
%   u       ~ (k*l) x 2         displacement field [u_x, u_y]
%   omega   ~ 1 x 4             evaluation region
%                                   [omega_1, omega_2] x [omega_3, omega_4]
%   s_grid  ~ 2 x 1             evaluation grid resolution = [k, l]
% OUT:
%   img_u   ~ k x l             interpolation of img over displaced grid
%   dimg_u  ~ (k*l) x (k*l*2)   gradient of interpol(img)
%--------------------------------------------------------------------------

% get image resolution
[m, n] = size(img);

% default grid resolution = image resolution
if nargin < 5
    s_grid = [m, n];
end
k = s_grid(1);      l = s_grid(2);

% default image region = [0, m * h_x] x [0, n * h_y]
if nargin < 4
    omega = [0, m * h(1), 0, n * h(2)];
end

% get cell centered grid over omega with resolution s_grid
[x, y] = cell_centered_grid(omega, [k, l]);
p = [x(:), y(:)];

% displace grid by u
p_displaced = p + u;

% interpolate img over resulting grid
[img_u, dimg_u] = bilinear_interpolation(img, h, p_displaced);
dimg_u = [spdiags(dimg_u(:, 1), 0, k*l, k*l), ...
    spdiags(dimg_u(:, 2), 0, k*l, k*l)];
img_u = reshape(img_u, [k, l]);

end