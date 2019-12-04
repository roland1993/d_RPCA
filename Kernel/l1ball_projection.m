function y = l1ball_projection(x)
%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------
% IN:
%   x ~ n x 1       input point
% OUT:
%   y ~ n x 1       projection of x onto {||.|| <= 1}
%--------------------------------------------------------------------------
% Linear Time Algorithm from
%  'Efficient Projections onto the l1-Ball for Learning in High Dimensions'
% by
%  John Duchi, Shai Shalev-Shwartz, Yoram Singer & Tushar Chandra
%--------------------------------------------------------------------------

x = x(:);       % force column vector
v = abs(x);     % precompute pointwise abs(x_i)
n = numel(x);

% catch easy case, where ||x||_1 <= 1 (no need for projection)
if sum(v) <= 1
    y = x;
else
    % reduction to simplex projection
    
    % initialization
    s = 0;
    rho = 0;
    U = true(n, 1);
    
    while any(U)
        
        % choose (random) index from U
        k = find(U, 1);
        
        % partition U
        G = (v >= v(k)) & U;
        L = (~G) & U;
        
        % compute s- and rho-updates
        d_rho = sum(G);
        d_s = sum(v(G));
        
        % update of U depending on s and rho
        if (s + d_s) - (rho + d_rho) * v(k) < 1
            s = s + d_s;
            rho = rho + d_rho;
            U = L;
        else
            U = G;
            U(k) = false;
        end
        
    end
    
    % compute final output
    theta = (s - 1) / rho;
    w = max(v - theta, 0);
    y = sign(x) .* w;
    
end

end