function e = matrix_norm(S)
%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------
% estimate matrix 2-norm
%--------------------------------------------------------------------------

tol = 1e-6;
maxIter = 100; % set max number of iterations to avoid infinite loop

i = 0;
e = zeros(maxIter, 1);

% get random starting point
x = randn(size(S, 2), 1);
x = x / norm(x);

while (i < 2) || (abs(e(i) - e(i - 1)) > tol * e(i))
    
    % increase iteration counter
    i = i + 1;
    
    % compute y = S'*S*x
    Sx = S*x;
    y = S'*Sx;
    
    % get || y ||_2
    normY = norm(y);
    
    % normalize y as input for next iteration
    x = y / normY;
    
    e(i) = sqrt(normY);
    
    % terminate at maxIter
    if i == maxIter, break; end
    
end

% delete redundant entries from e
e(i + 1 : maxIter) = [];

end