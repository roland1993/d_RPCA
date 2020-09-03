%   min_L delta_{|| . || <= nu}(B * L)
%       + sum_i || l_i - I_(u_i) ||_1
%       + sum_i mu * CURVATURE(u_i)
%       + delta_{mean(u_x) = 0, mean(u_y) = 0}
%
%   MEAN-FREE & NO REFERENCE & USES UNIQUENESS-TERM

function [u0, L0, SV_history] = mf_nn_curvature_registration_no_ref_ml(img, optPara)
%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------
% IN:
%   img     ~ cell(k, 1)        array of images
%   optPara ~ struct            optimization parameters with fields
%       .theta      ~ 1 x 1     over-relaxation parameter
%       .maxIter    ~ 1 x 1     maximum number of iterations
%       .tol        ~ 1 x 1     tolerance for p/d-gap + infeasibilities
%       .outerIter  ~ 1 x 1     number of outer iterations
%       .mu         ~ 1 x 1     weighting factor (see model above)
%       .nu_factor  ~ 1 x 1     reduction of nu per outer iterate
%       .bc         ~ string    boundary condition for grid discretization
%       .doPlots    ~ logical   do plots during optimization?
%
% OUT:
%   u0              ~ m*n x 2 x k        	displacement fields
%   L0              ~ m*n x k   			low-rank components
%	SV_history 		~ cell					development of singular values
%--------------------------------------------------------------------------

% make sure that interpolation routines are on search path
if ~exist('evaluate_displacement.m', 'file')
    addpath(genpath('..'));
end

% some local function handles
vec = @(x) x(:);
normalize = @(x) (x - min(x(:))) / (max(x(:)) - min(x(:)));

% get number of template images
k = length(img);

% normalize images to range [0, 1]
for i = 1 : k
    img{i} = normalize(img{i});
end

% get image resolution etc.
[m, n] = size(img{1});
omega = [0, m, 0, n];

%-------BEGIN MULTI LEVEL-------------------------------------------------%
% resolution at lowest level: m, n >= 2 ^ 5
numLevels = min(floor(log2([m, n]) - 5)) + 1;

% initialize output
% u = cell(numLevels, max(optPara.outerIter));
% L = cell(numLevels, max(optPara.outerIter));
SV_history = cell(numLevels, 1);

% get multi-level representation
ML = cell(numLevels, k);
for i = 1 : k
    
    % highest level = input
    ML{numLevels, i} = img{i};
    
    % downsampling
    for lev = (numLevels - 1) : (-1) : 1
        ML{lev, i} = conv2(ML{lev + 1, i}, 0.25 * ones(2), 'same');
        ML{lev, i} = ML{lev, i}(1 : 2 : end, 1 : 2 : end);
    end
    
end

% iterate over levels
for lev = 1 : numLevels
    
    % get image resolution at current level
    [m, n] = size(ML{lev, 1});
    
    % image region
    h_grid = (omega([2, 4]) - omega([1, 3])) ./ [m, n];
    
    % set optimization parameters
    theta = optPara.theta;
    tol = optPara.tol;
    mu = optPara.mu;
    bc = optPara.bc;
    doPlots = optPara.doPlots;
    maxIter = optPara.maxIter;
    if lev == 1
        outerIter = optPara.outerIter(1);
        nu_factor = optPara.nu_factor(1);
    else
        outerIter = optPara.outerIter(2);
        nu_factor = optPara.nu_factor(2);
    end
    
    % initialize primal and dual variables (x and p)
    if lev == 1
        x = zeros(3 * k * m * n, 1);
    else
        
        % get resolution from last level
        [m_low_res, n_low_res] = size(ML{lev - 1, 1});
        
        % prolongation of primal variables to higher resolution
        x_high_res = zeros(2 * m_low_res, 2 * n_low_res, 3 * k);
        x_low_res = reshape(x, m_low_res, n_low_res, 3 * k);
        for i = 1 : (3 * k)
            x_high_res(:, :, i) = kron(x_low_res(:, :, i), ones(2));
        end
        x = vec(x_high_res(1 : m, 1 : n, :));
        
    end
    p = zeros(4 * k * m * n, 1);
    
    % gradient operator for displacement fields
    A2 = discrete_laplacian(m, n, h_grid, k, bc);
    
    % all zeros
    A3 = sparse(k * m * n, 2 * k * m * n);
    
    % identity matrix
    A4 = speye(k * m * n);
    
    % all zeros
    A5 = sparse(2 * k * m * n, k * m * n);
    
    % mean free operator
    A6 = mean_free_operator(m, n, k);
    
    % set function handle to G-part of target function
    G_handle = @(x, c_flag) G(x, [m, n, k], c_flag);
    
    % if plotting was requested -> create figure (handle)
    if doPlots && (lev == 1)
        fh1 = figure;
        fh2 = figure;
        fh3 = figure;
    end
    
    % initialize displacement field u0
    u0 = reshape(x(1 : 2 * k * m * n), m * n, 2, k);
    
    % initialize storage for singular values
    SV_history{lev} = zeros(k, max(outerIter));

    %-------BEGIN OUTER ITERATION-----------------------------------------%
    for o = 1 : outerIter
        
        % reference vector for computing SAD from L
        b = zeros(k * m * n, 1);
        
        % templates evaluated using current u0
        T_u = zeros(m, n, k);
        dT = cell(k, 1);
        for i = 1 : k
            [T_u(:, :, i), dT{i}] = evaluate_displacement( ...
                ML{lev, i}, h_grid, u0(:, :, i), omega);
            b((i - 1) * m * n + 1 : i * m * n) = ...
                vec(T_u(:, :, i)) - dT{i} * vec(u0(:, :, i));
        end
        
        % on first iterate: estimate nu from nn mean-free image-matrix
        if (o == 1) && (lev == 1)
            D = reshape(A6 * vec(T_u), m * n, k);
            [~, S, ~] = svd(D, 'econ');
            nu = nu_factor * sum(diag(S));
        elseif (o == 1) && (lev > 1)
            nu = 2 * nu_factor * nu;
        else
            nu = nu_factor * nu;
        end
        
        % upper left block of A ~> image gradients
        A1 = -blkdiag(dT{:});
        
        % build up A from the computed blocks
        A = [       A1,     A4
                    A2,     A5
                    A3,     A6          ];
        
        % estimate spectral norm of A
        e = matrix_norm(A);
        norm_A_est = e(end);
        
        % use estimated norm to get primal and dual step sizes
        tau = sqrt(0.99 / norm_A_est ^ 2);
        sigma = sqrt(0.99 / norm_A_est ^ 2);
        
        % get function handle to F-part of target function
        F_handle = @(y, c_flag) F(y, b, k, h_grid, mu, nu, sigma, c_flag);
        
        % perform optimization
        [x, p, primal_history, dual_history] = chambolle_pock( ...
            F_handle, G_handle, A, x, p, theta, tau, sigma, maxIter, tol);
        
        % get displacements and low rank components from minimizer x
        u0 = reshape(x(1 : 2 * k * m * n), m * n, 2, k);
        L0 = x(2 * k * m * n + 1 : end);
        
        % store info on development of singular values
        [~, S, ~] = svd(reshape(A6 * L0, m * n, k), 'econ');
        SV_history{lev}(:, o) = diag(S);
        
        % plot progress (if requested)
        if doPlots
            plot_progress(fh1, primal_history, dual_history);
            set(fh1, 'NumberTitle', 'off', ...
                'Name', sprintf('ITERATE %d OUT OF %d', o, outerIter));
            display_results(ML(lev, :), u0, [], reshape(L0, m, n, k), fh2);
            plot_sv(fh3, SV_history{lev}(:, 1 : o));
            drawnow;
        end

    end
    %-------END OUTER ITERATION-------------------------------------------%
    
end
%-------END MULTI LEVEL---------------------------------------------------%


%-------BEGIN LOCAL FUNCTION DEFINITIONS----------------------------------%
    function [res1, res2, res3] = ...
            F(y, b, k, h_grid, mu, nu, sigma, conjugate_flag)
        % splits input y = [y1; y2; y3] and computes
        %   F_1(y1) = h1 * h2 * || y1 - b ||_1
        %   F_2(y2) = 0.5 * mu * h1 * h2 * || y2 ||_2^2
        %   F_3(y3) = delta_{|| . ||_* <= nu}(y3)
        
        % get number of template images and number of pixels per image
        mn = numel(y) / (4 * k);
        
        % split input y into r- and v-part
        y1 = y(1 : k * mn);
        y2 = y(k * mn + 1 : 3 * k * mn);
        y3 = y((3 * k * mn) + 1 : end);
        
        if nargout == 3
            
            % initialize output
            res3 = zeros(4 * k * mn, 1);
            
            % apply SAD to y1-part
            [~, ~, res3_F1] = ...
                SAD(y1, b, prod(h_grid), sigma, conjugate_flag);
            res3(1 : k * mn) = res3_F1;
            
            % apply sum of squares to y2
            g = sparse(2 * k * mn, 1);
            [~, ~, res3_F2] = ...
                SSD(y2, g, mu * prod(h_grid), sigma, conjugate_flag);
            res3(k * mn + 1 : 3 * k * mn) = res3_F2;
            
            % apply delta_{|| . ||_* <= nu} to y3
            [~, ~, res3_F3] = nuclear_norm_constraint( ...
                y3, k, sigma, nu, conjugate_flag);
            res3(3 * k * mn + 1 : end) = res3_F3;
            
            % dummy outputs
            res1 = [];
            res2 = [];
            
        else
            
            % apply F1 = SAD to y1-part
            [res1_F1, res2_F1] = ...
                SAD(y1, b, prod(h_grid), sigma, conjugate_flag);
            
            % apply sum of squares to y2
            g = sparse(2 * k * mn, 1);
            [res1_F2, res2_F2] = ...
                SSD(y2, g, mu * prod(h_grid), sigma, conjugate_flag);
            
            % apply delta_{|| . ||_* <= nu} to y3
            [res1_F3, res2_F3] = nuclear_norm_constraint( ...
                y3, k, sigma, nu, conjugate_flag);
            
            % compute outputs
            res1 = [res1_F1, res1_F2, res1_F3];
            res2 = max([res2_F1, res2_F2, res2_F3]);
            
        end
        
    end
%-------------------------------------------------------------------------%
    function [res1, res2, res3] = G(x, s, conjugate_flag)
        
        % s = [m, n, k]
        
        % split x = [x_u, x_l]
        x_u = x(1 : 2 * prod(s));
        x_l = x(2 * prod(s) + 1 : end);
        
        if nargout == 3
            
            % apply delta_{mean(u_x) = 0, mean(u_y) = 0} to x_u
            [~, ~, res3_G1] = ...
                mean_zero_indicator(x_u, s, conjugate_flag);
            
            % apply zero-function to x_l
            [~, ~, res3_G2] = zero_function(x_l, conjugate_flag);
            
            % combine outputs
            res3 = [res3_G1; res3_G2];
            
            % dummy outputs
            res1 = [];
            res2 = [];
            
        else
            
            % apply delta_{mean(u_x) = 0, mean(u_y) = 0} to x_u
            [res1_G1, res2_G1] = ...
                mean_zero_indicator(x_u, s, conjugate_flag);
            
            % apply zero-function to x_l
            [res1_G2, res2_G2] = zero_function(x_l, conjugate_flag);
            
            % combine outputs
            res1 = res1_G1 + res1_G2;
            res2 = max([res2_G1, res2_G2]);
            
        end
        
    end
%-------------------------------------------------------------------------%
    function plot_progress(fh, primal_history, dual_history)
        % plots progress of one outer iterate
        
        % make figure active
        figure(fh);

        % plot primal vs. dual energy
        subplot(2, 2, 1);
        plot(primal_history(:, 1));
        hold on;
        plot(dual_history(:, 1));
        hold off;
        axis tight;
        grid on;
        xlabel('#iter');
        legend({'primal energy', 'dual energy'}, ...
            'Location', 'SouthOutside', ...
            'Orientation', 'Horizontal');
        title('primal vs. dual')
        
        % plot numerical gap
        GAP = abs((primal_history(:, 1) - dual_history(:, 1)) ./ ...
            dual_history(:, 1));
        subplot(2, 2, 2);
        semilogy(GAP);
        axis tight;
        grid on;
        xlabel('#iter');
        legend({'absolute primal-dual gap'}, ...
            'Location', 'SouthOutside', ...
            'Orientation', 'Horizontal');
        title('primal-dual gap');
        
        % plot constraints
        subplot(2, 2, 3);
        semilogy(primal_history(:, 6));
        hold on;
        semilogy(primal_history(:, 7));
        semilogy(dual_history(:, 6));
        semilogy(dual_history(:, 7));
        hold off;
        axis tight;
        grid on;
        xlabel('#iter');
        legend({'F', 'G', 'F*', 'G*'}, ...
            'Location', 'SouthOutside', 'Orientation', 'Horizontal');
        title('constraints');
        
        % plot different components of F
        subplot(2, 2, 4);
        plot(primal_history(:, 1));
        hold on;
        plot(primal_history(:, 2));
        plot(primal_history(:, 3));
        hold off;
        axis tight;
        ylim([0, max(primal_history(:, 1))]);
        grid on;
        xlabel('#iter');
        legend({'F', '\Sigma_i || T_i(u_i) - l_i ||_1', ...
            '\Sigma_i TV(u_i)'}, 'Location', 'SouthOutside', ...
            'Orientation', 'Horizontal');
        title('decomposition of F');
        
    end
%-------------------------------------------------------------------------%
    function plot_sv(fh, SV_history)
        
        % get number of images = number of sigular values
        numImg = size(SV_history, 1);
        
        % get numImg colors
        cmap = jet(numImg);
        
        % cell-array for legend entries
        names = cell(numImg + 1, 1);
        
        % do the plotting
        figure(fh);
        clf;
        hold on;
        for j = 1 : numImg
            plot(SV_history(j, :), '-x', 'Color', cmap(j, :));
            names{j} = ['\sigma_{', num2str(j), '}'];
        end
        plot(sum(SV_history, 1), 'k--x');
        names{numImg + 1} = '\Sigma_i \sigma_i';
        hold off;
        xlim([0.5, size(SV_history, 2) + 0.5]);
        xlabel('#outer iter');
        title('development of singular values');
        grid on;
        legend(names);
        
    end
%-------END LOCAL FUNCTION DEFINITIONS------------------------------------%

end