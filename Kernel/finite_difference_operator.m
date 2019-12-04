function D = finite_difference_operator(m, n, h_grid, k, bc)
%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------
% IN:
%   m       ~ 1 x 1                         number of rows
%   n       ~ 1 x 1                         number of columns
%   h_grid  ~ 1 x 2                         grid spacing
%   k       ~ 1 x 1                         number of grids
%   bc      ~ string                        boundary condition
% OUT:
%   D       ~ sparse(4*k*m*n x 2*k*m*n)     finite difference operator for
%                                               k displacement fields in 
%                                               x- and y-direction
%--------------------------------------------------------------------------

% some standard parameters
if nargin < 5, bc = 'linear'; end
if nargin < 4, k = 1; end

% x-derivative
Dx = (1 / h_grid(1)) * spdiags([-ones(m, 1), ones(m, 1)], 0 : 1, m, m);
if strcmp(bc, 'linear')
    Dx(m, [m - 1, m]) = [-1, 1];
elseif strcmp(bc, 'neumann')
    Dx(m, m) = 0;
else
    error('Unknown boundary condition!');
end

% y-derivative
Dy = (1 / h_grid(2)) * spdiags([-ones(n, 1), ones(n, 1)], 0 : 1, n, n);
if strcmp(bc, 'linear')
    Dy(n, [n - 1, n]) = [-1, 1];
elseif strcmp(bc, 'neumann')
    Dy(n, n) = 0;
else
    error('Unknown boundary condition!');
end

% construct operator for one grid and extend to k grids
G = kron(speye(2), [kron(speye(n), Dx); kron(Dy, speye(m))]);
D = kron(speye(k), G);

end