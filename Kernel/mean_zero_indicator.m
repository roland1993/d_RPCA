function [res1, res2, res3] = mean_zero_indicator(u, s, conjugate_flag)
%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------
% IN:
%       u               ~ 2*m*n*k x 1   vector of k displacement fields
%       s               ~ 3 x 1         dimensions & number of images
%       conjugate_flag  ~ logical       evalutate indicator or conjugate?
% OUT:
%   IF ~conjugate_flag
%       res1            ~ 1 x 1         delta_{mean([u_x; u_y]) = 0}(u)
%       res2            ~ 1 x 1         constraint violation measure
%       res3            ~ 2*m*n*k x 1   prox = projection onto <v, r> = 0
%   IF conjugate_flag
%       res1            ~ 1 x 1         delta_{span{r}}(u)
%       res2            ~ 1 x 1         constraint violation measure
%       res3            ~ 2*m*n*k x 1   prox = projection onto span{r}
%--------------------------------------------------------------------------

% fetch images dimensions etc.
m = s(1);
n = s(2);
k = s(3);

% seperate x- and y-components
x_idx = repmat([true(m * n, 1); false(m * n, 1)], [k, 1]);
u_x = u(x_idx);
y_idx = repmat([false(m * n, 1); true(m * n, 1)], [k, 1]);
u_y = u(y_idx);

% normal-vector of mean-zero subspace
r = ones(m * n * k, 1);
norm_r_squared = m * n * k;

if ~conjugate_flag
    
    if nargout == 3
        
        % initialize prox
        res3 = zeros(size(u));
        
        % projection of u_x to subspace r' * v = 0    <=>   u_x-mean = 0
        res3(x_idx) = u_x - ((r' * u_x) / norm_r_squared) * r;
        
        % projection of u_y to subspace r' * v = 0    <=>   u_y-mean = 0
        res3(y_idx) = u_y - ((r' * u_y) / norm_r_squared) * r;
        
        % dummy outputs
        res1 = [];
        res2 = [];
        
    else
        
        % fctn. value for indicator
        res1 = 0;
        
        % distance from mean = 0
        res2 = max(abs([mean(u_x), mean(u_y)]));
        
    end
    
else
    
    % fctn. value for indicator
    res1 = 0;
    
    % initialize prox
    res3 = zeros(size(u));
    
    % projection of u_x to subspace span{r}
    res3(x_idx) = ((r' * u_x) / norm_r_squared) * r;
    
    % projection of u_y to subspace span{r}
    res3(y_idx) = ((r' * u_x) / norm_r_squared) * r;
    
    % constraint measure: distance to span{r}
    res2 = max([abs(u_x - res3(x_idx)); abs(u_y - res3(y_idx))]);
    
end

end