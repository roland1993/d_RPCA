function [img_p, dimg_p] = bilinear_interpolation(img, h, p)
%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------
% IN:
%   img         ~ m x n     image
%   h           ~ 2 x 1     image step sizes [h_x, h_y]
%   p           ~ l x 2     evaluation points (in xy-coordinates)
% OUT:
%   img_p       ~ l x 1     bilinear interpolation of img at points p
%   dimg_p      ~ l x 2     x-/y-derivative of interpol(img) at points p
%--------------------------------------------------------------------------

% get some more parameters
[m, n] = size(img);
l = size(p, 1);

% use homogeneous coordinates for p
p = [p, ones(size(p(:, 1)))];

% define world matrix to switch from xy-coordinates to ij-coordinates
W = [1/h(1), 0, 1/2; ...
    0, 1/h(2), 1/2; ...
    0, 0, 1];

% switch coordinates (and remove homogeneous component)
q = W * p';
q = q(1 : 2, :)';

% evaluate bilinear interpolation (and x/y-derivative) of img at points p
img_p = zeros(l, 1);
dimg_p = zeros(l, 2);

% get top left, bottom left, top right and bottom right grid point
sub_tl = floor(q);
sub_bl = sub_tl + [1, 0];
sub_tr = sub_tl + [0, 1];
sub_br = sub_tl + [1, 1];

% get valid subscripts (i.e. inside region [1, m] x [1, n])
tl_valid = (1 <= sub_tl(:, 1)) & (sub_tl(:, 1) <= m) & ...
    (1 <= sub_tl(:, 2)) & (sub_tl(:, 2) <= n);
bl_valid = (1 <= sub_bl(:, 1)) & (sub_bl(:, 1) <= m) & ...
    (1 <= sub_bl(:, 2)) & (sub_bl(:, 2) <= n);
tr_valid = (1 <= sub_tr(:, 1)) & (sub_tr(:, 1) <= m) & ...
    (1 <= sub_tr(:, 2)) & (sub_tr(:, 2) <= n);
br_valid = (1 <= sub_br(:, 1)) & (sub_br(:, 1) <= m) & ...
    (1 <= sub_br(:, 2)) & (sub_br(:, 2) <= n);

% convert all valid subscripts to linear indices
idx_tl = sub2ind([m, n], sub_tl(tl_valid, 1), sub_tl(tl_valid, 2));
idx_bl = sub2ind([m, n], sub_bl(bl_valid, 1), sub_bl(bl_valid, 2));
idx_tr = sub2ind([m, n], sub_tr(tr_valid, 1), sub_tr(tr_valid, 2));
idx_br = sub2ind([m, n], sub_br(br_valid, 1), sub_br(br_valid, 2));

% get weighting coefficients
chi = q - sub_tl;

% work corners piece by piece
%   1. top left
img_p(tl_valid) = img_p(tl_valid) + ...
    (1 - chi(tl_valid, 1)) .* (1 - chi(tl_valid, 2)) .* img(idx_tl);

dimg_p(tl_valid, 1) = dimg_p(tl_valid, 1) + (1 / h(1)) * ...
    (-1) * (1 - chi(tl_valid, 2)) .* img(idx_tl);

dimg_p(tl_valid, 2) = dimg_p(tl_valid, 2) + (1 / h(2)) * ...
    (-1) * (1 - chi(tl_valid, 1)) .* img(idx_tl);

%   2. bottom left
img_p(bl_valid) = img_p(bl_valid) + ...
    chi(bl_valid, 1) .* (1 - chi(bl_valid, 2)) .* img(idx_bl);

dimg_p(bl_valid, 1) = dimg_p(bl_valid, 1) + (1 / h(1)) * ...
    (+1) * (1 - chi(bl_valid, 2)) .* img(idx_bl);

dimg_p(bl_valid, 2) = dimg_p(bl_valid, 2) + (1 / h(2)) * ...
    (-1) * chi(bl_valid, 1) .* img(idx_bl);

%   3. top right
img_p(tr_valid) = img_p(tr_valid) + ...
    (1 - chi(tr_valid, 1)) .* chi(tr_valid, 2) .* img(idx_tr);

dimg_p(tr_valid, 1) = dimg_p(tr_valid, 1) + (1 / h(1)) * ...
    (-1) * chi(tr_valid, 2) .* img(idx_tr);

dimg_p(tr_valid, 2) = dimg_p(tr_valid, 2) + (1 / h(2)) * ...
    (+1) * (1 - chi(tr_valid, 1)) .* img(idx_tr);

%   4. bottom right
img_p(br_valid) = img_p(br_valid) + ...
    chi(br_valid, 1) .* chi(br_valid, 2) .* img(idx_br);

dimg_p(br_valid, 1) = dimg_p(br_valid, 1) + (1 / h(1)) * ...
    (+1) * chi(br_valid, 2) .* img(idx_br);

dimg_p(br_valid, 2) = dimg_p(br_valid, 2) + (1 / h(2)) * ...
    (+1) * chi(br_valid, 1) .* img(idx_br);

end