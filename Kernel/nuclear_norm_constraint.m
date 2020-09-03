function [res1, res2, res3] = ...
    nuclear_norm_constraint(L, numImg, tau, nu, conjugate_flag)
%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------
% IN:
%       L               ~ m*n*numImg x 1    all images in one column vector
%       numImg          ~ 1 x 1             number of images
%       tau             ~ 1 x 1             prox step size
%       nu              ~ 1 x 1             constraint threshold
%       conjugate_flag  ~ logical           eval. constraint or conjugate?
% OUT:
%   IF NOT conjugate_flag:
%       res1            ~ 1 x 1             indicator of ||.||_* <= nu
%       res2            ~ 1 x 1             constraint violation measure
%       res3            ~ m*n x 1           prox-step of indicator
%   IF conjugate_flag:
%       res1            ~ 1 x 1             weighted spectral norm
%       res2            ~ 1 x 1             constraint violation measure
%       res3            ~ m*n x 1           prox-step of spectral norm
%--------------------------------------------------------------------------

% by default: evaluate constraint instead of its conjugate
if nargin < 5, conjugate_flag = false; end

% reshape L back into a matrix ~ m*n x numImg
L = reshape(L, [], numImg);

% compute svd of L
[U, S, V] = svd(L, 'econ');
S = diag(S);

if ~conjugate_flag
    
    % get prox operator (l1-ball-projection of SV-vector)
    if nargout == 3
        
        res3 = U * diag(nu * l1ball_projection(S / nu)) * V';
        res3 = res3(:);
        
        % dummy outputs
        res1 = [];
        res2 = [];
        
    else
        
        % distance to feasible region
        if sum(S) > nu
            res2 = (sum(S) - nu) / nu;
        else
            res2 = 0;
        end
        
        % return 0 as fctn. value of indicator
        res1 = 0;
        
    end
    
else
    
    % compute prox of spectral norm via prox of inf-norm of sv-vector S
    if nargout == 3
        
        % nu and tau in one factor
        mu = nu * tau;
        
        % prox on S via Moreau's identity
        %   -> conjugate of inf-norm is l1-ball indicator
        S_prox = S - mu * l1ball_projection(S / mu);
        
        % use prox step of sv-vector to compute prox of spectral norm
        res3 = U * diag(S_prox) * V';
        res3 = res3(:);
        
        % dummy outputs
        res1 = [];
        res2 = [];
        
    else
        
        % conjugate of nn-constraint = spectral norm = max singular value
        res1 = nu * max(S);
        
        % no constraints present
        res2 = 0;
        
    end
    
end

end