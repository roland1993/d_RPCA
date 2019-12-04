function [res1, res2, res3] = SAD_weight(L, I, mu, sigma, conjugate_flag)
%--------------------------------------------------------------------------
% This file is part of my master's thesis entitled
%           'Low rank- and sparsity-based image registration'
% For the whole project see
%           https://github.com/roland1993/MA
% If you have questions contact me at
%           roland.haase [at] student.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------
% IN:
%       L               ~ m*n*numImg x 1    input variables
%       I               ~ m*n*numImg x 1    I from SAD(L) = ||L - I||_1
%       mu              ~ 1 x 1             weighting factor
%       sigma           ~ 1 x 1             prox step size
%       conjugate_flag  ~ logical           evaluate SAD or conjugate?
% OUT:
%   IF NOT conjugate_flag:
%       res1            ~ 1 x 1             function value SAD(L)
%       res2            ~ 1 x 1             constraint violation measure
%       res3            ~ m*n x 1           prox-step of SAD for L
%   IF conjugate_flag:
%       res1            ~ 1 x 1             convex conjugate SAD*(L)
%       res2            ~ 1 x 1             constraint violation measure
%       res3            ~ m*n x 1           prox-step of SAD* for L
%--------------------------------------------------------------------------

% by default: evaluate SAD instead of its conjugate
if nargin < 5, conjugate_flag = false; end

if ~conjugate_flag
    % EITHER ~> evaluate SAD and Prox_[SAD] at L
    
    if nargout == 3
        
        % dummy outputs
        res1 = [];
        res2 = [];
        
        % prox-step for ||L - I||_1 =: SAD ~> pointwise shrinkage
        res3 = zeros(size(L));
        diff_LI = L - I;
        idx1 = (diff_LI > (sigma * mu));
        idx2 = (diff_LI < (-1) * (sigma * mu));
        idx3 = ~(idx1 | idx2);
        res3(idx1) = L(idx1) - (sigma * mu);
        res3(idx2) = L(idx2) + (sigma * mu);
        res3(idx3) = I(idx3);
        
    else
        
        % compute SAD
        res1 = mu * sum(abs(L - I));
        
        % constraint measure;
        res2 = 0;
        
    end
    
else
    % OR ~> evaluate SAD* and Prox_[SAD*] at L
    
    % compute prox-step for SAD* with Moreau's identity (if requested)
    if nargout == 3
        
        % dummy outputs
        res1 = [];
        res2 = [];
        
        [~, ~, prox] = SAD_weight(L / sigma, I, mu, 1 / sigma, false);
        res3 = L - sigma * prox;
        
    else
        
        % compute SAD*(L) = delta_{||.||_inf <= 1}(L) + <L,I>
        if max(abs(L)) > mu
            res2 = max(abs(L)) - mu;
        else
            res2 = 0;
        end
        
        res1 = mu * ((L / mu)' * I);
        
    end
    
end

end