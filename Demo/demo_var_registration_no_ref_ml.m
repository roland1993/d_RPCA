%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------

% demo script for var_registration_no_ref_ml.m

% add project to path
addpath(genpath('..'));

% clean-up
clear all, close all, clc;

% generate data
m = 200;    n = 200;    k = 10;
[data, LM_data] = dynamicTestImage(m, n, k);
img = cell(k, 1);
LM = cell(k, 1);
for i = 1 : k
    img{i} = data(:, :, i);
    LM{i} = LM_data(:, :, i);
end
omega = [0, m, 0, n];

% set optimization parameters
optPara.theta = 1;
optPara.maxIter = 2000;
optPara.tol = 1e-3;
optPara.outerIter = [16 2];
optPara.mu = 1e-1;
optPara.bc = 'neumann';
optPara.doPlots = true;

% call registration routine
u = var_registration_no_ref_ml(img, optPara);

% fetch results
uStar = u{end, optPara.outerIter(2)};

% evaluate results
img_u = cell(k, 1);
LM_transformed = cell(k, 1);

for i = 1 : k
    
    % get transformed images
    img_u{i} = evaluate_displacement(img{i}, [1, 1], uStar(:, :, i));
    
    % transform landmarks
    LM_transformed{i} = ...
        landmark_transform(LM{i}, reshape(uStar(:, :, i), [m n 2]), omega);
    
end

% landmark accuracy in terms of mean distance to mean lm position
LM_acc = landmark_accuracy(LM);
LM_transformed_acc = landmark_accuracy(LM_transformed);

%% display results

% input, output and low rank components in comparison
figure;
colormap gray(256);

while true
    for i = 1 : k
        
        subplot(1, 2, 1);
        imshow(img{i}, [], 'InitialMagnification', 'fit');
        hold on
        scatter(LM{i}(:, 2), LM{i}(:, 1));
        hold off
        title(sprintf('input T_{%d}', i));
        
        subplot(1, 2, 2);
        imshow(img_u{i}, [], 'InitialMagnification', 'fit');
        hold on
        scatter(LM_transformed{i}(:, 2), LM_transformed{i}(:, 1));
        hold off
        title(sprintf('output T_{%d}(u_{%d})', i, i));
        
        waitforbuttonpress;
        
    end
end