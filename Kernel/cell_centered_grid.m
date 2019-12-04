function [x, y] = cell_centered_grid(omega, s)
%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------
% IN:
%   omega   ~ 1 x 4     grid region [omega_1, omega_2] x [omega_3, omega_4]
%   s       ~ 2 x 1     grid size [m, n]
% OUT:
%   x       ~ m x n     x-coordinates of grid points
%   y       ~ m x n     y-coordinates of grid points
%--------------------------------------------------------------------------

m = s(1);   n = s(2);
h_x = (omega(2) - omega(1)) / m;
h_y = (omega(4) - omega(3)) / n;

x = repmat(...
    omega(1) + (h_x * ((1/2) : 1 : (m - (1/2))))', ...
    [1, n]);
y = repmat(...
    omega(3) + (h_y * ((1/2) : 1 : (n - (1/2)))), ...
    [m, 1]);

end