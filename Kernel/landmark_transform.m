function LM_transformed = landmark_transform(LM, u, omega)
%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------
% IN:
%   LM              ~ numLM x 2     landmarks in original sytem
%   u               ~ m x n x 2     deformation field
%   omega           ~ 1 x 4         image region
% OUT:
%   LM_transformed  ~ numLM x 2     landmarks in transformed sytem
%--------------------------------------------------------------------------

% get number of landmarks
numLM = size(LM, 1);

% get resolution of u
m = size(u, 1);
n = size(u, 2);

% get grid widths
hx = (omega(2) - omega(1)) / m;
hy = (omega(4) - omega(3)) / n;

%
F = @(z) z + [bilinear_interpolation(u(:, :, 1), [hx, hy], z), ...
    bilinear_interpolation(u(:, :, 2), [hx, hy], z)];

% initialize output
LM_transformed = 0 * LM;

% generate cell centered grid over omega
[xx, yy] = cell_centered_grid(omega, [m, n]);
p = [xx(:), yy(:)];

% transformed grid
g = p + reshape(u, [m * n, 2]);

for i = 1 : numLM
    
    % find closest point in transformed grid to current landmark
    %   -> initial guess for transformed lm
    [~, min_idx]  = min(sum((g - LM(i, :)) .^ 2, 2));
    
    % fixed-point iteration for landmark inversion
    y = LM(i, :);
    x = p(min_idx, :);
    for j = 1 : 100
        x = x + (y - F(x));
        if sum(F(x) - LM(i, :) .^ 2) < 1e-12
            break;
        end
    end
    if sum(F(x) - LM(i, :) .^ 2) >= 1e-12
        warning('Fixed-point iteration failed! Returning initial guess.');
        x = p(min_idx, :);
    end
    LM_transformed(i, :) = x;
    
end

end