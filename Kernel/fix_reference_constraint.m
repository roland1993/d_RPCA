function [res1, res2, res3] = ...
    fix_reference_constraint(u, s, ref_idx, conjugate_flag)
%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------
% IN:
%   u                   ~ 2*k*m*n x 1   k displacement fields
%   s                   ~ 3 x 1         s = [m, n, k]
%   ref_idx             ~ 1 x 1         index of reference image
%   conjugate_flag      ~ string        evaluate constraint or prox?
% OUT:
%   IF conjugate_flag:
%       res1            ~ 1 x 1         delta_{||.||_{u^ref_idx=0}
%       res2            ~ 1 x 1         constraint violation measure
%       res3            ~ m*n*4 x 1     prox-step
%   IF NOT conjugate_flag:
%       res1            ~ 1 x 1         conjugate evaluated at u
%       res2            ~ 1 x 1         constraint violation measure
%       res3            ~ m*n*4 x 1     prox-step
%--------------------------------------------------------------------------


m = s(1);
n = s(2);
k = s(3);

u = reshape(u, 2 * m * n, k);

res1 = zeros(1);
res2 = zeros(1);
res3 = zeros(size(u));

for i = 1 : k
    
    if i == ref_idx
        [res1_i, res2_i, res3_i] = zero_function(u(:, i), ~conjugate_flag);
    else
        [res1_i, res2_i, res3_i] = zero_function(u(:, i), conjugate_flag);
    end
    
    res1 = res1 + res1_i;
    res2 = max(res2, res2_i);
    res3(:, i) = res3_i;
    
end

res3 = res3(:);

end