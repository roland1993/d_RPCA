function [res1, res2, res3] = norm21(v, mu, sigma, conjugate_flag)
%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------
% IN:
%       v               ~ m*n*4 x 1     input vector
%       mu              ~ 1 x 1         weighting parameter
%       sigma           ~ 1 x 1         prox step size
%       conjugate_flag  ~ logical       eval. norm or its conjugate?
% OUT:
%   IF conjugate_flag:
%       res1            ~ 1 x 1         delta_{||.||_{2,inf} <= mu}(v)
%       res2            ~ 1 x 1         constraint violation measure
%       res3            ~ m*n*4 x 1     prox of indicator delta_{...}
%   IF NOT conjugate_flag:
%       res1            ~ 1 x 1         mu * ||v||_{2,1}
%       res2            ~ 1 x 1         constraint violation measure
%       res3            ~ m*n*4 x 1     prox_[mu * ||.||_{2,1}](v)
%--------------------------------------------------------------------------

% by default: evaluate ||v||_{2,1} instead of its conjugate
if nargin < 4, conjugate_flag = false; end

if ~conjugate_flag
    % EITHER ~> evaluate mu*||.||_{2,1} and its prox at v
    
    if nargout == 3
        
        % compute prox-step for F(.) = ||.||_{2,1} with Moreau's identity
        %   [(id + sigma * dF)^(-1)](v) =
        %       v - sigma * [(id + (1 / sigma) * dF*)^(-1)](v / sigma)
        [~, ~, conj_prox] = norm21(v(:) / sigma, mu, 1 / sigma, true);
        res3 = v(:) - sigma * conj_prox;
        
        % dummy outputs
        res1 = [];
        res2 = [];
        
    else
        
        % constraint measure = 0
        res2 = 0;
        
        % reshape v into 4 columns and compute pointwise 2-norm
        v = reshape(v, [], 4);
        norm_v = sqrt(sum(v .^ 2, 2));
        
        % mu * ||v||_{2,1} = mu * sum_i ||v_i||_2
        res1 = mu * sum(norm_v);
        
    end
    
else
    % OR ~> evaluate delta_{||.||_{2,inf} <= mu} and its prox at v
    
    % reshape v into 4 columns and compute pointwise 2-norm
    v = reshape(v, [], 4);
    norm_v = sqrt(sum(v .^ 2, 2));
    
    if nargout == 3
        
        % pointwise prox: v_i := (mu * v_i) / max(mu, ||v_i||_2)
        n = max(norm_v, mu);
        res3 = (mu * v) ./ n;
        res3 = res3(:);
        
        % dummy outputs
        res1 = [];
        res2 = [];
        
    else
        
        % conjugate [mu * ||.||_{2,1}]* = delta_{||.||_{2,inf} <= mu}
        if max(norm_v) > mu
            res2 = max(norm_v) - mu;
        else
            res2 = 0;
        end
        
        % fctn. value of indicator
        res1 = 0;

    end
    
end

end