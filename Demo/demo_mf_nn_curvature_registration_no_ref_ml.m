%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------

% demo script for mf_nn_curvature_registration_no_ref_ml.m
clear all, close all, clc;

%% choose dataset from {heart, kidney}
dataset = 'kidney';

switch dataset
    
    case 'heart'
        
        % load data
        load('heart_mri.mat');
        k = length(IDX);
        
        % downsampling
        img = cell(k, 1);
        factor = 2;
        for i = 1 : k
            img{i} = conv2(data(:, :, IDX(i)), ...
                ones(factor) / factor ^ 2, 'same');
            img{i} = img{i}(1 : factor : end, 1 : factor : end);
        end
        [m, n] = size(img{1});
        omega = [0, m, 0, n];
        
        % adjust scaling of landmarks
        LM = cell(1, k);
        for i = 1 : k
            LM{i} = [m n] .* LM_IDX{i};
        end
        
        % set optimization parameters
        optPara.theta = 1;
        optPara.maxIter = 2000;
        optPara.tol = 5e-3;
        optPara.outerIter = [16 2];
        optPara.mu = 1e1;
        optPara.bc = 'neumann';
        optPara.nu_factor = [0.95 0.95];
        optPara.doPlots = false;
        
    case 'kidney'
        
        % load data
        load('dcemri_kidney.mat');
        
        k = size(data, 3);
        m = size(data, 1);
        n = size(data, 2);
        omega = [0, m, 0, n];
        
        img = cell(k, 1);
        for i = 1 : k
            img{i} = data(:, :, i);
            img{i} = (img{i} - min(img{i}(:))) / ...
                (max(img{i}(:)) - min(img{i}(:)));
            LM{i} = [m n] .* LM{i};
        end
        
        % set optimization parameters
        optPara.theta = 1;
        optPara.maxIter = 2000;
        optPara.tol = 5e-3;
        optPara.outerIter = [16 2];
        optPara.mu = 2e1;
        optPara.bc = 'dirichlet';
        optPara.nu_factor = [0.95 0.95];
        optPara.doPlots = false;
        
    otherwise
        error('No such dataset!');
        
end

%% call registration routine
tic;
[uStar, L, SV_history] = ...
    mf_nn_curvature_registration_no_ref_ml(img, optPara);
runtime_registration = toc;

%% evaluate results
LStar = reshape(L, m, n, k);
img_u = cell(k, 1);

LM_available = exist('LM', 'var');
if LM_available
    LM_transformed = cell(k, 1);
end
    
for i = 1 : k
    
    % get transformed images
    img_u{i} = evaluate_displacement(img{i}, [1, 1], uStar(:, :, i));
    
    % transform landmarks
    if LM_available
        LM_transformed{i} = landmark_transform(LM{i}, ...
            reshape(uStar(:, :, i), [m n 2]), omega);
    end
        
end

% landmark accuracy in terms of mean distance to mean lm position
if LM_available
    LM_acc = landmark_accuracy(LM);
    LM_transformed_acc = landmark_accuracy(LM_transformed);
end

%% save results
method = mfilename;
if ~exist('results', 'file')
    mkdir('.', 'results');
end
save(sprintf('./results/res_%s_mf_nn_curvature', dataset));

%% display results
if optPara.doPlots
    
    [cc_x, cc_y] = cell_centered_grid(omega, [m n]);
    yStar = reshape(uStar, [m n 2 k]) + ...
        repmat(cat(3, cc_x, cc_y), [1 1 1 k]);
    green = cat(3, zeros(m, n), ones(m, n), zeros(m, n));
    
    figure;
    colormap gray(256);
    
    while true
        for i = 1 : k
            
            subplot(1, 3, 1);
            imagesc(img{i});
            axis image;
            hold on;
            if LM_available
                scatter(LM{i}(:, 2), LM{i}(:, 1), ...
                    'bo', 'MarkerFaceColor', 'red');
            end
            plot(yStar(1 : 3 : end, 1 : 3 : end, 2, i), ...
                yStar(1 : 3 : end, 1 : 3 : end, 1, i), 'g-');
            plot(yStar(1 : 3 : end, 1 : 3 : end, 2, i)', ...
                yStar(1 : 3 : end, 1 : 3 : end, 1, i)', 'g-');
            hold off;
            title(sprintf('input T_{%d}', i));
            
            subplot(1, 3, 2);
            imagesc(img_u{i});
            axis image;
            if LM_available
                hold on
                scatter(LM_transformed{i}(:, 2), ...
                    LM_transformed{i}(:, 1), ...
                    'bo', 'MarkerFaceColor', 'red');
                hold off
                for j = 1 : size(LM_transformed{i}, 1)
                    text(LM_transformed{i}(j, 2), ...
                        LM_transformed{i}(j, 1), sprintf('  %d', j), ...
                        'Color', 'green');
                end
            end
            title(sprintf('output T_{%d}(u_{%d})', i, i));
            
            subplot(1, 3, 3);
            imshow(LStar(:, :, i) - mean(LStar, 3), [], ...
                'InitialMagnification', 'fit');
            hold on;
            imagesc(...
                'YData', omega(1) + [0.5, m - 0.5], ...
                'XData', omega(3) + [0.5, n - 0.5], ...
                'CData', green, ...
                'AlphaData', abs(img_u{i} - LStar(:, :, i)));
            hold off;
            title(sprintf('output L_{%d}', i));
            
            waitforbuttonpress;
            
        end
    end
    
end