function [res1, res2, res3] = SSD(x, g, lambda, tau, conjugate_flag)
%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------
% IN:
%       x               ~ m*n*numImg x 1    input variables
%       g               ~ m*n*numImg x 1    g in SSD(x) = 0.5||x - g||_2^2
%       lambda          ~ 1 x 1             weighting parameter
%       tau             ~ 1 x 1             prox step size
%       conjugate_flag  ~ logical           evaluate conjugate?
% OUT:
%   IF NOT conjugate_flag:
%       res1            ~ 1 x 1             function value lambda*SSD(x)
%       res2            ~ 1 x 1             constraint violation measure
%       res3            ~ m*n x 1           prox of lambda*SSD for x
%   IF conjugate_flag:
%       res1            ~ 1 x 1             conjugate [lambda*SSD]*(x)
%       res2            ~ 1 x 1             constraint violation measure
%       res3            ~ m*n x 1           prox of [lambda*SSD]* for x
%--------------------------------------------------------------------------

% % by default: evaluate SSD instead of its conjugate
% if nargin < 4, conjugate_flag = false; end

% initialize constraint violation measure with 0
res2 = 0;

if ~conjugate_flag
    % EITHER ~> evaluate SSD and Prox_[SSD] at x
    
    if nargout == 3
        
        res1 = [];
        
        % prox-step
        res3 = (x + (lambda * tau) * g) / (1 + (lambda * tau));

    else
        
        % compute SSD
        res1 = lambda * 0.5 * sum((x - g) .^ 2);
        
    end
    
else
    % OR ~> evaluate SSD* and Prox_[SSD*] at x
    
    % compute prox-step for SSD* with Moreau's identity (if requested)
    if nargout == 3
        
        res1 = [];
        
        [~, ~, prox] = SSD(x / tau, g, lambda, 1 / tau, false);
        res3 = x - tau * prox;
        
    else
        
        % compute SSD*(x) = 0.5 * ||x||_2^2 + <x, g>
        res1 = ...
            lambda * (0.5 * sum((x / lambda) .^ 2) + ((x / lambda)' * g));
        
    end
    
end

end