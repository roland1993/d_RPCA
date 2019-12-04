function img_u = display_results(img, u, refIdx, L, fh)
%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------

%-----------------INITIALIZATION------------------------------------------%
% get image resolution etc.
[m, n] = size(img{1});
h_img = [1, 1];
omega = [0, m, 0, n];
h_grid = (omega([2, 4]) - omega([1, 3])) ./ [m, n];

% get cc-grid for plotting
[cc_x, cc_y] = cell_centered_grid(omega, [m, n]);
cc_grid = [cc_x(:), cc_y(:)];

% get green image for visualizing image differences
green = cat(3, zeros(m, n), ones(m, n), zeros(m, n));

% do the plotting
if exist('fh', 'var') && ~isempty(fh)
    figure(fh);
    colormap gray(256);
else
    fh = figure;
    colormap gray(256);
end

set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);

%-----------------EVALUATE DISPLACEMENTS----------------------------------%

% set flag for models with/without reference
reference = exist('refIdx', 'var') && ~isempty(refIdx);

% set flag for models with/without low-rank-part
lowrank = exist('L', 'var') && ~isempty(L);

if reference
    
    % k = number of templates
    k = length(img) - 1;
    
    % get reference image
    R = img{refIdx};
    
    % some indexing
    tempIdx = 1 : (k + 1);
    tempIdx(refIdx) = [];
    IDX = [tempIdx, refIdx];
    
    % evaluate displacements
    img_u = cell(k + 1, 1);
    for i = 1 : k
        img_u{tempIdx(i)} = ...
            evaluate_displacement(img{tempIdx(i)}, h_img, u(:, :, i));
    end
    img_u{refIdx} = img{refIdx};
    
else
    
    % k = number of templates
    k = length(img);
    
    % evaluate displacements
    img_u = cell(k, 1);
    for i = 1 : k
        img_u{i} = ...
            evaluate_displacement(img{i}, h_img, u(:, :, i));
    end
    
end

% indexing for selected quivers
if m >= 50
    IDX1 = round(linspace(1, m, 50));
else
    IDX1 = 1 : m;
end
if n >= 50
    IDX2 = round(linspace(1, n, 50));
else
    IDX2 = 1 : n;
end
IDX = (IDX2 - 1) * m + IDX1';
IDX = IDX(:);

%-----------------PLOTTING------------------------------------------------%

% case 1: with reference & low-rank-part
if reference && lowrank
    
    % get mean of L
    meanL = mean(L, 3);
    
    for i = 1 : (k + 1)
        
        subplot(3, k + 1, i);
        imshow(img{IDX(i)}, [0 1], 'InitialMagnification', 'fit');
        if i <= k
            hold on;
            quiver(cc_grid(IDX, 2), cc_grid(IDX, 1), ...
                u(IDX, 2, i), u(IDX, 1, i), 1, 'r');
            hold off;
            title(sprintf('T_{%d} with u_{%d}', i, i));
        else
            title('R');
        end
        
        subplot(3, k + 1, (k + 1) + i);
        imshow(img_u{IDX(i)}, [0 1], 'InitialMagnification', 'fit');
        hold on;
        imagesc(...
            'YData', omega(1) + h_grid(1) * [0.5, m - 0.5], ...
            'XData', omega(3) + h_grid(2) * [0.5, n - 0.5], ...
            'CData', green, ...
            'AlphaData', abs(img_u{IDX(i)} - L(:, :, i)));
        hold off;
        if i <= k
            title(sprintf( ...
                'T_{%d}(u_{%d}) with |T_{%d}(u_{%d}) - l_{%d}|', ...
                i, i, i, i, i));
        else
            title(sprintf('R with |R - l_{%d}|', i));
        end
        
        subplot(3, k + 1, 2 * (k + 1) + i);
        imshow(L(:, :, i) - meanL, [], 'InitialMagnification', 'fit');
        title(sprintf('l_{%d} - l_{mean}', i));
        
    end
    
end
    
% case 2: with reference & without low-rank-part
if reference && ~lowrank
    
    for i = 1 : (k + 1)
        
        subplot(2, k + 1, i);
        imshow(img{IDX(i)}, [0 1], 'InitialMagnification', 'fit');
        if i <= k
            hold on;
            quiver(cc_grid(IDX, 2), cc_grid(IDX, 1), ...
                u(IDX, 2, i), u(IDX, 1, i), 1, 'r');
            hold off;
            title(sprintf('T_{%d} with u_{%d}', i, i));
        else
            title('R');
        end
        
        if i <= k
            subplot(2, k + 1, (k + 1) + i);
            imshow(img_u{IDX(i)}, [0 1], 'InitialMagnification', 'fit');
            hold on;
            imagesc(...
                'YData', omega(1) + h_grid(1) * [0.5, m - 0.5], ...
                'XData', omega(3) + h_grid(2) * [0.5, n - 0.5], ...
                'CData', green, ...
                'AlphaData', abs(img_u{IDX(i)} - R));
            hold off;
            title(sprintf('T_{%d}(u_{%d}) with |T_{%d}(u_{%d}) - R|', ...
                i, i, i, i));
        end
        
    end
    
end

% case 3: without reference & with low-rank-part
if ~reference && lowrank
    
    % get mean of L
    meanL = mean(L, 3);
    
    for i = 1 : k
        
        subplot(3, k, i);
        imshow(img{i}, [0 1], 'InitialMagnification', 'fit');
        hold on;
        quiver(cc_grid(IDX, 2), cc_grid(IDX, 1), ...
            u(IDX, 2, i), u(IDX, 1, i), 1, 'r');
        hold off;
        title(sprintf('T_{%d} with u_{%d}', i, i));
        
        subplot(3, k, k + i);
        imshow(img_u{i}, [0 1], 'InitialMagnification', 'fit');
        hold on;
        imagesc(...
            'YData', omega(1) + h_grid(1) * [0.5, m - 0.5], ...
            'XData', omega(3) + h_grid(2) * [0.5, n - 0.5], ...
            'CData', green, ...
            'AlphaData', abs(img_u{i} - L(:, :, i)));
        hold off;
        title(sprintf('T_{%d}(u_{%d}) with |T_{%d}(u_{%d}) - l_{%d}|', ...
            i, i, i, i, i));
        
        subplot(3, k, 2 * k + i);
        imshow(L(:, :, i) - meanL, [], 'InitialMagnification', 'fit');
        title(sprintf('l_{%d} - l_{mean}', i));
        
    end
    
end

% case 4: without reference & without low-rank-part
if ~reference && ~lowrank
    
    % get mean image for plotting difference
    meanImg = zeros(size(img_u{1}));
    for i = 1 : k
        meanImg = meanImg + img_u{i} / k;
    end
    
    for i = 1 : k
        
        subplot(2, k, i);
        imshow(img{i}, [0 1], 'InitialMagnification', 'fit');
        hold on;
        quiver(cc_grid(IDX, 2), cc_grid(IDX, 1), ...
            u(IDX, 2, i), u(IDX, 1, i), 1, 'r');
        hold off;
        title(sprintf('T_{%d} with u_{%d}', i, i));
        
        subplot(2, k, k + i);
        imshow(img_u{i}, [0 1], 'InitialMagnification', 'fit');
        hold on;
        imagesc(...
            'YData', omega(1) + h_grid(1) * [0.5, m - 0.5], ...
            'XData', omega(3) + h_grid(2) * [0.5, n - 0.5], ...
            'CData', green, ...
            'AlphaData', abs(img_u{i} - meanImg));
        hold off;
        title(sprintf( ...
            'T_{%d}(u_{%d}) with |T_{%d}(u_{%d}) - T_{mean}|', ...
            i, i, i, i));
        
    end
    
end

end