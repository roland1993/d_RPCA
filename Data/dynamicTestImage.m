function [T, LM] = dynamicTestImage(m, n, numFrames)
%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------

if ~exist('m', 'var'), m = 200; end
if ~exist('n', 'var'), n = 200; end
if ~exist('numFrames', 'var'), numFrames = 6; end

[xx, yy] = meshgrid(linspace(-1, 1, n), linspace(-1, 1, m));
T = zeros(m, n, numFrames);
LM = zeros(16, 2, numFrames);

f1 = 4 * pi;
p1 = 0.375 * pi;
f2 = 6 * pi;
p2 = 0.125 * pi;

for i = 1 : numFrames
    
    dx = -0.1 * sin(pi * (i) / numFrames);
    dy = 0.1 * cos(pi * (i) / numFrames);
    
    ellipse_rad = 0.4;
    ellipse = double( ...
        sqrt(2 * (xx + dx - 0.25) .^ 2 + (yy + dy) .^ 2) <= ellipse_rad);
    if mod(i, 2) == 0
        texture = sin(f1 * (yy + dy) + p1) .^ 2;
    else
        texture = sin(f2 * (xx + dx) + p2) .^ 2;
    end
    IDX = sqrt(2 * (xx + dx - 0.25) .^ 2 + (yy + dy) .^ 2) ...
            <= 0.6 * ellipse_rad;
    ellipse(IDX) = texture(IDX);
    
    frame_rad = 0.7;
    frame_width = 0.15;
    tmp = 1;
    frame = double(reshape( ...
            frame_rad <= max(abs(tmp * [xx(:), yy(:)]), [], 2) & ...
            max(abs(tmp * [xx(:), yy(:)]), [], 2) <= ...
                                frame_rad + frame_width, m, n));
    
    rect = double(reshape((-0.6 <= xx(:)) & (xx(:) <= -0.2) & ...
        (-0.5 <= yy(:)) & (yy(:) <= 0.5), m, n));

    tmp = imgaussfilt(rect + ellipse + frame, (m + n) / 150);
    T(:, :, i) = tmp;
    
    % set landmarks
    LM(1, :, i) = [-0.5, -0.6];
    LM(2, :, i) = [-0.5, -0.2];
    LM(3, :, i) = [0.5, -0.6];
    LM(4, :, i) = [0.5, -0.2];
    LM(5, :, i) = [-dy, (0.4 / sqrt(2)) - (dx - 0.25)];
    LM(6, :, i) = [-dy, -(0.4 / sqrt(2)) - (dx - 0.25)];
    LM(7, :, i) = [(0.4 - dy), (0.25 - dx)];
    LM(8, : ,i) = [(-0.4 - dy), (0.25 - dx)];
    LM(9, :, i) = [(-dy), (0.25 - dx)];
    LM(10, :, i) = [0.7, 0.7];
    LM(11, :, i) = [0.7, -0.7];
    LM(12, :, i) = [-0.7, 0.7];
    LM(13, :, i) = [-0.7, -0.7];
    LM(14, :, i) = [0.85, 0.85];
    LM(15, :, i) = [0.85, -0.85];
    LM(16, :, i) = [-0.85, 0.85];
    LM(17, :, i) = [-0.85, -0.85];
    
end

% convert landmarks to coordinate system omega = [0, m] x [0, n]
LM = LM + 1;
LM(:, 1, :) = LM(:, 1, :) * (m / 2);
LM(:, 2, :) = LM(:, 2, :) * (n / 2);

% figure;
% colormap gray(256);
% for i = 1 : numFrames
%     
%     imagesc(T(:, :, i));
%     axis image;
%     colorbar;
%     drawnow;
%     pause(1 / numFrames);
%     
% end

end