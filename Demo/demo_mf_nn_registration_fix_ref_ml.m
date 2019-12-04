%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------

% demo script for mf_nn_registration_fix_ref_ml.m
% INSTRUCTIONS:
%   1. Request CDSOF-datasets from 
%       https://hci.iwr.uni-heidelberg.de/benchmarks/Challenging_Data_for_Stereo_and_Optical_Flow
%   2. Unpack into Data-directory (i.e., '../Data').
%   3. Get sub_imread.m from
%       https://hci.iwr.uni-heidelberg.de/system/files/private/datasets/256921177/sub_imread.zip
%      and add it to a path known to MATLAB.

% add project to path
addpath(genpath('..'));

% clean-up
clear all, close all, clc;

% normalization to [0, 1]
normalize = @(y) (y - min(y(:))) / (max(y(:)) - min(y(:)));

% select data_set from {'BlinkingArrow', 'FlyingSnow', 'ShadowOnTruck'}
data_set = 'ShadowOnTruck';
data_path = sprintf('../Data/ChallengingSequences/%s/sequence', data_set);

switch data_set
    
    case 'FlyingSnow'
        
        % load + downsample images
        IDX = 1 : 10;
        k = numel(IDX);
        img = cell(1, k);
        factor = 2;
        for i = 1 : k
            img{i} = normalize(double( ...
                sub_imread(sprintf('%s/%06d_0.pgm', data_path, IDX(i)))));
            img{i} = conv2(img{i}, ones(factor) / factor ^ 2, 'same');
            img{i} = img{i}(1 : factor : end, 1 : factor : end);
        end
        
        % optimization parameters
        optPara.theta = 1;
        optPara.maxIter = 2000;
        optPara.tol = 1e-3;
        optPara.outerIter = [16 2];
        optPara.mu = 4.5e-2;
        optPara.nu_factor = [0.91 0.91];
        optPara.bc = 'neumann';
        optPara.doPlots = true;
        
        % select reference
        ref_idx = 4 - IDX(1) + 1;
        
    case 'ShadowOnTruck'
        
        % load + downsample images
        IDX = 0 : 30;
        k = numel(IDX);
        img = cell(1, k);
        factor = 2;
        for i = 1 : k
            img{i} = normalize(double( ...
                sub_imread(sprintf('%s/%06d_0.pgm', data_path, IDX(i)))));
            img{i} = conv2(img{i}, ones(factor) / factor ^ 2, 'same');
            img{i} = img{i}(1 : factor : end, 1 : factor : end);
        end
        
        % optimization parameters
        optPara.theta = 1;
        optPara.maxIter = 2000;
        optPara.tol = 1e-3;
        optPara.outerIter = [16 2];
        optPara.mu = 4.5e-2;
        optPara.nu_factor = [0.95 0.95];
        optPara.bc = 'neumann';
        optPara.doPlots = true;
        
        % select reference
        ref_idx = 13 - IDX(1) + 1;
        
    case 'BlinkingArrow'
        
        % load + downsample images
        IDX = 1 : 10;
        k = numel(IDX);
        img = cell(1, k);
        factor = 2;
        for i = 1 : k
            img{i} = normalize(double( ...
                sub_imread(sprintf('%s/%06d_0.pgm', data_path, IDX(i)))));
            img{i} = conv2(img{i}, ones(factor) / factor ^ 2, 'same');
            img{i} = img{i}(1 : factor : end, 1 : factor : end);
        end
        
        % optimization parameters
        optPara.theta = 1;
        optPara.maxIter = 2000;
        optPara.tol = 1e-3;
        optPara.outerIter = [16 2];
        optPara.mu = 7.5e-2;
        optPara.nu_factor = [0.91 0.91];
        optPara.bc = 'neumann';
        optPara.doPlots = true;
        
        % select reference
        ref_idx = 3 - IDX(1) + 1;
        
end

% call registration routine
[u, L, SV_history] = mf_nn_registration_fix_ref_ml(img, ref_idx, optPara);

% evaluate results
uStar = u{end, optPara.outerIter(2)};
LStar = L{end, optPara.outerIter(2)};
img_u = cell(k, 1);
for i = 1 : k
    img_u{i} = evaluate_displacement(img{i}, [1, 1], uStar(:, :, i));
end

%% display results

% input, output and low rank components in comparison
figure;
colormap gray(256);
while true
    for i = 1 : k
        
        subplot(1, 3, 1);
        imshow(img{i}, [], 'InitialMagnification', 'fit');
        title(sprintf('input T_{%d}', i));
        
        subplot(1, 3, 2);
        imshow(img_u{i}, [], 'InitialMagnification', 'fit');
        title(sprintf('output T_{%d}(u_{%d})', i, i));
        
        subplot(1, 3, 3);
        imshow(LStar(:, :, i), [], 'InitialMagnification', 'fit');
        title(sprintf('output L_{%d}', i));
        
        waitforbuttonpress;
        
    end
end